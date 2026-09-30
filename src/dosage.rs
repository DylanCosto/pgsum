//! Quantitative diploid genotypes. Decimal text is preserved without binary floating-point rounding.
use crate::decimal::Decimal;
use crate::pack::Weight;
use crate::score::{Coefficient, Contribution, ExactSum, contribution};
use crate::term::Model;
use num_bigint::BigInt;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "UPPERCASE")]
pub enum Field {
    #[default]
    Gt,
    Ds,
    Gp,
}

impl Field {
    pub fn label(self) -> &'static str {
        match self {
            Self::Gt => "GT",
            Self::Ds => "DS",
            Self::Gp => "GP",
        }
    }
}

/// ALT dosage, or probabilities for 0, 1 and 2 copies of one ALT allele.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "field", content = "value")]
pub enum Measurement {
    DS(Decimal),
    GP([Decimal; 3]),
}

fn amount(d: &Decimal) -> Contribution {
    Contribution {
        negative: d.negative,
        coefficient: Coefficient::Big(d.coefficient.parse().expect("validated decimal")),
        exponent: d.exponent,
    }
}

fn bounded(text: &str, max: u8) -> Option<Decimal> {
    let d = Decimal::parse(text)?;
    // Compare exact integers at a common scale; f64 can round an out-of-range value into range.
    let exponent = d.exponent.min(0);
    let magnitude = d.coefficient.parse::<BigInt>().ok()? * BigInt::from(10u8).pow((d.exponent - exponent) as u32);
    let limit = BigInt::from(max) * BigInt::from(10u8).pow((-exponent) as u32);
    ((!d.negative || d.coefficient == "0") && magnitude <= limit).then_some(d)
}

fn sum(values: impl IntoIterator<Item = Contribution>) -> Decimal {
    let mut total = ExactSum::default();
    for value in values {
        total.add(&value);
    }
    total.finish()
}

fn product(c: &Contribution, d: &Decimal) -> Contribution {
    let left = match &c.coefficient {
        Coefficient::Small(v) => BigInt::from(*v),
        Coefficient::Big(v) => v.clone(),
    };
    Contribution {
        negative: c.negative ^ d.negative,
        coefficient: Coefficient::Big(left * d.coefficient.parse::<BigInt>().expect("validated decimal")),
        exponent: c.exponent + d.exponent,
    }
}

impl Measurement {
    pub fn parse(field: Field, text: &str) -> Option<Self> {
        match field {
            Field::Gt => None,
            Field::Ds => Some(Self::DS(bounded(text, 2)?)),
            Field::Gp => {
                let values: Vec<_> = text.split(',').map(|s| bounded(s, 1)).collect::<Option<_>>()?;
                let probabilities: [Decimal; 3] = values.try_into().ok()?;
                // Explicit strict policy: no silent renormalization or clipping of probabilities.
                let total = sum(probabilities.iter().map(amount)).rescaled(0);
                (total.coefficient == "1" && total.exponent == 0).then_some(Self::GP(probabilities))
            }
        }
    }

    /// Revalidate serialized measurements, including public Decimal fields.
    pub fn valid(&self) -> bool {
        let safe = |d: &Decimal| {
            d.coefficient.len() <= crate::decimal::MAX_TEXT
                && !d.coefficient.is_empty()
                && d.coefficient.bytes().all(|c| c.is_ascii_digit())
                && (-1127..=100).contains(&d.exponent)
                && (d.coefficient == "0" || !d.coefficient.starts_with('0'))
        };
        if !match self {
            Self::DS(d) => safe(d),
            Self::GP(p) => p.iter().all(safe),
        } {
            return false;
        }
        match self {
            Self::DS(d) => Self::parse(Field::Ds, &d.to_python_string()).as_ref() == Some(self),
            Self::GP(p) => {
                Self::parse(
                    Field::Gp,
                    &p.iter().map(Decimal::to_python_string).collect::<Vec<_>>().join(","),
                )
                .as_ref()
                    == Some(self)
            }
        }
    }

    pub fn alt_dosage(&self) -> Decimal {
        match self {
            Self::DS(d) => d.clone(),
            Self::GP(p) => {
                let mut twice = amount(&p[2]);
                twice.coefficient = Coefficient::Big(p[2].coefficient.parse::<BigInt>().expect("decimal") * 2);
                sum([amount(&p[1]), twice])
            }
        }
    }

    pub fn effect_dosage(&self, effect_is_alt: bool) -> Decimal {
        let dosage = self.alt_dosage();
        if effect_is_alt {
            dosage
        } else {
            sum([amount(&Decimal::parse("2").unwrap()), amount(&dosage).negated()])
        }
    }

    /// DS identifies an additive expectation only. GP supports expectations of all supported models.
    pub fn contribution(&self, model: Model, weights: &[Weight], effect_is_alt: bool) -> Option<Contribution> {
        match self {
            Self::DS(_) if model == Model::Additive => Some(product(
                &contribution(model, weights, 1)?,
                &self.effect_dosage(effect_is_alt),
            )),
            Self::DS(_) => None,
            Self::GP(p) => {
                let values = (0..3)
                    .map(|i| {
                        let dosage = if effect_is_alt { i } else { 2 - i };
                        Some(product(&contribution(model, weights, dosage as u8)?, &p[i]))
                    })
                    .collect::<Option<Vec<_>>>()?;
                Some(amount(&sum(values)))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn score(field: Field, input: &str, model: Model, weights: &[&str], alt: bool) -> String {
        Measurement::parse(field, input)
            .unwrap()
            .contribution(
                model,
                &weights.iter().map(|s| Weight::Text(s.to_string())).collect::<Vec<_>>(),
                alt,
            )
            .unwrap()
            .to_decimal()
            .to_plain_string()
    }
    #[test]
    fn decimal_dosages_and_reversed_alleles_are_exact() {
        assert_eq!(
            score(Field::Ds, "0.1234567890123456789", Model::Additive, &["-0.2"], true),
            "-0.02469135780246913578"
        );
        assert_eq!(score(Field::Ds, "0.125", Model::Additive, &["0.2"], false), "0.3750");
        assert!(
            Measurement::parse(Field::Ds, "0.5")
                .unwrap()
                .contribution(Model::Dominant, &[Weight::Text("1".into())], true)
                .is_none()
        );
    }
    #[test]
    fn probabilities_use_expected_model_contributions() {
        assert_eq!(score(Field::Gp, "0.2,0.3,0.5", Model::Additive, &["2"], true), "2.6");
        assert_eq!(score(Field::Gp, "0.2,0.3,0.5", Model::Dominant, &["2"], true), "1.6");
        assert_eq!(score(Field::Gp, "0.2,0.3,0.5", Model::Recessive, &["2"], false), "0.4");
        assert_eq!(
            score(Field::Gp, "0.2,0.3,0.5", Model::DosageWeights, &["-1", "2", "4"], true),
            "2.4"
        );
    }
    #[test]
    fn rejects_invalid_range_and_probability_mass() {
        for s in ["-0.1", "2.0000000000000000001", "NaN", "inf", ".", "0.5,0.2"] {
            assert!(Measurement::parse(Field::Ds, s).is_none(), "{s}");
        }
        for s in [
            "0.2,0.3",
            "0.2,0.3,0.4",
            "0.2,0.3,0.5000000000000000001",
            "-0.2,0.7,0.5",
            "0,0,0",
        ] {
            assert!(Measurement::parse(Field::Gp, s).is_none(), "{s}");
        }
    }
    #[test]
    fn extreme_accepted_exponents_do_not_underflow() {
        let value = score(
            Field::Ds,
            "1.0000000000e-1000",
            Model::Additive,
            &["1.0000000000e-1000"],
            true,
        );
        assert_eq!(value, format!("0.{}1{}", "0".repeat(1999), "0".repeat(20)));
    }
}
