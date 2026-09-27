//! Exact decimal weights.
//!
//! Scoring files give weights as decimal text. They are parsed exactly (no floating point) into a coefficient
//! and a power-of-ten exponent, following the same acceptance rules and the same text form as Python's
//! `decimal.Decimal`, so results can be compared as strings.

/// A finite decimal: `(-1)^negative × coefficient × 10^exponent`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decimal {
    pub negative: bool,
    /// ASCII digits without leading zeros (`"0"` for zero).
    pub coefficient: String,
    pub exponent: i64,
}

/// Longest accepted weight text.
pub const MAX_TEXT: usize = 128;
/// Accepted range of the adjusted exponent (the exponent of the leading digit).
pub const ADJUSTED_RANGE: std::ops::RangeInclusive<i64> = -1000..=100;

impl Decimal {
    /// Parse `[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?` of at most 128 bytes whose adjusted exponent is in
    /// −1000..=100; anything else is `None`.
    pub fn parse(text: &str) -> Option<Decimal> {
        if text.len() > MAX_TEXT {
            return None;
        }
        let bytes = text.as_bytes();
        let mut i = 0;
        let negative = match bytes.first() {
            Some(b'-') => {
                i += 1;
                true
            }
            Some(b'+') => {
                i += 1;
                false
            }
            _ => false,
        };
        let digits = |from: usize| from + bytes[from..].iter().take_while(|b| b.is_ascii_digit()).count();
        let int_end = digits(i);
        let (frac_start, frac_end) = if bytes.get(int_end) == Some(&b'.') {
            (int_end + 1, digits(int_end + 1))
        } else {
            (int_end, int_end)
        };
        if int_end == i && frac_end == frac_start {
            return None;
        }
        let mut exponent: i64 = 0;
        let mut end = frac_end;
        if matches!(bytes.get(end), Some(b'e' | b'E')) {
            let mut j = end + 1;
            let exp_negative = match bytes.get(j) {
                Some(b'-') => {
                    j += 1;
                    true
                }
                Some(b'+') => {
                    j += 1;
                    false
                }
                _ => false,
            };
            let exp_end = digits(j);
            if exp_end == j {
                return None;
            }
            // An exponent too large for i64 is far outside the accepted range anyway.
            let magnitude: i64 = text[j..exp_end].parse().ok()?;
            exponent = if exp_negative { -magnitude } else { magnitude };
            end = exp_end;
        }
        if end != bytes.len() {
            return None;
        }
        let all_digits = format!("{}{}", &text[i..int_end], &text[frac_start..frac_end]);
        let trimmed = all_digits.trim_start_matches('0');
        let coefficient = if trimmed.is_empty() {
            "0".to_owned()
        } else {
            trimmed.to_owned()
        };
        let exponent = exponent.checked_sub((frac_end - frac_start) as i64)?;
        let decimal = Decimal {
            negative,
            coefficient,
            exponent,
        };
        ADJUSTED_RANGE.contains(&decimal.adjusted()).then_some(decimal)
    }

    /// Exponent of the leading digit.
    /// The same value written with `exponent`: trailing zeros removed or appended. Removing digits that are
    /// not zeros would change the value, so they are kept (the exponent then stays finer).
    pub fn rescaled(&self, exponent: i64) -> Decimal {
        if self.coefficient == "0" {
            return Decimal {
                negative: false,
                coefficient: "0".into(),
                exponent,
            };
        }
        let mut coefficient = self.coefficient.clone();
        let mut e = self.exponent;
        while e < exponent && coefficient.ends_with('0') && coefficient.len() > 1 {
            coefficient.pop();
            e += 1;
        }
        while e > exponent {
            coefficient.push('0');
            e -= 1;
        }
        Decimal {
            negative: self.negative,
            coefficient,
            exponent: e,
        }
    }

    pub fn adjusted(&self) -> i64 {
        self.exponent + self.coefficient.len() as i64 - 1
    }

    /// The coefficient as a signed 64-bit integer, when it fits.
    pub fn mantissa_i64(&self) -> Option<i64> {
        let magnitude: i64 = self.coefficient.parse().ok()?;
        Some(if self.negative { -magnitude } else { magnitude })
    }

    /// Python's `format(Decimal, "f")` form: never an exponent.
    pub fn to_plain_string(&self) -> String {
        let c = &self.coefficient;
        let sign = if self.negative { "-" } else { "" };
        if self.exponent >= 0 {
            let zeros = if c == "0" { 0 } else { self.exponent as usize };
            return format!("{sign}{c}{}", "0".repeat(zeros));
        }
        let places = (-self.exponent) as usize;
        let padded = format!("{}{c}", "0".repeat((places + 1).saturating_sub(c.len())));
        let (int_part, frac_part) = padded.split_at(padded.len() - places);
        format!("{sign}{int_part}.{frac_part}")
    }

    /// Python's `str(Decimal)` form.
    pub fn to_python_string(&self) -> String {
        let c = &self.coefficient;
        let left_digits = self.exponent + c.len() as i64;
        let dot_place = if self.exponent <= 0 && left_digits > -6 {
            left_digits
        } else {
            1
        };
        let (int_part, frac_part) = if dot_place <= 0 {
            ("0".to_owned(), format!(".{}{}", "0".repeat((-dot_place) as usize), c))
        } else if dot_place as usize >= c.len() {
            (
                format!("{}{}", c, "0".repeat(dot_place as usize - c.len())),
                String::new(),
            )
        } else {
            (
                c[..dot_place as usize].to_owned(),
                format!(".{}", &c[dot_place as usize..]),
            )
        };
        let exp = if left_digits == dot_place {
            String::new()
        } else {
            format!("E{:+}", left_digits - dot_place)
        };
        format!(
            "{}{}{}{}",
            if self.negative { "-" } else { "" },
            int_part,
            frac_part,
            exp
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rescaling_keeps_the_value() {
        let d = |t: &str| Decimal::parse(t).unwrap();
        assert_eq!(d("1.2000").rescaled(-1).to_python_string(), "1.2");
        assert_eq!(d("1.2").rescaled(-4).to_python_string(), "1.2000");
        assert_eq!(d("0E-23").rescaled(0).to_python_string(), "0");
        // Digits that are not zeros stay.
        assert_eq!(d("1.25").rescaled(-1).to_python_string(), "1.25");
    }

    fn s(text: &str) -> String {
        Decimal::parse(text).unwrap().to_python_string()
    }

    #[test]
    fn accepts_and_formats_like_python() {
        // Expected values are Python's str(Decimal(text)).
        assert_eq!(s("0.16220387987485377"), "0.16220387987485377");
        assert_eq!(s("0012.30"), "12.30");
        assert_eq!(s("-0"), "-0");
        assert_eq!(s("0.00"), "0.00");
        assert_eq!(s("1e2"), "1E+2");
        assert_eq!(s("1.5E-7"), "1.5E-7");
        assert_eq!(s("0.000001"), "0.000001");
        assert_eq!(s("0.0000001"), "1E-7");
        assert_eq!(s(".5"), "0.5");
        assert_eq!(s("5."), "5");
        assert_eq!(s("+3"), "3");
        assert_eq!(s("123E-2"), "1.23");
        assert_eq!(s("-4.50e1"), "-45.0");
    }

    #[test]
    fn rejects_non_numbers_and_out_of_range() {
        for text in [
            "", ".", "e5", "1e", "1.2.3", "NaN", "inf", " 1", "1 ", "1_0", "0x10", "--1",
        ] {
            assert_eq!(Decimal::parse(text), None, "{text:?}");
        }
        assert!(Decimal::parse("1e100").is_some());
        assert!(Decimal::parse("1e101").is_none());
        assert!(Decimal::parse("1e-1000").is_some());
        assert!(Decimal::parse("1e-1001").is_none());
        assert!(Decimal::parse("0e-1001").is_none());
        assert!(Decimal::parse("1e99999999999999999999").is_none());
        assert!(Decimal::parse(&"1".repeat(129)).is_none());
    }

    #[test]
    fn mantissa() {
        let d = Decimal::parse("-0.16220387987485377").unwrap();
        assert_eq!((d.mantissa_i64(), d.exponent), (Some(-16220387987485377), -17));
    }
}
