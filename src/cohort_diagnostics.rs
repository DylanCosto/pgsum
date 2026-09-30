//! Compact missing-call states and source-record ordinals, loaded separately from packed cohort calls.
use crate::genotype::State;
use serde::{Deserialize, Serialize};

pub const PASSING: u8 = 255;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) enum MissingStates {
    #[default]
    None,
    Uniform(u8),
    Runs(Vec<(u32, u32, u8)>),
    Dense(Vec<u8>),
}

impl MissingStates {
    pub fn from_codes(states: Vec<u8>) -> Self {
        if states.iter().all(|s| *s == PASSING) {
            return Self::None;
        }
        if states.iter().all(|s| *s == states[0]) {
            return Self::Uniform(states[0]);
        }
        let mut runs = Vec::new();
        let mut i = 0;
        while i < states.len() {
            if states[i] == PASSING {
                i += 1;
                continue;
            }
            let start = i;
            while i < states.len() && states[i] == states[start] {
                i += 1;
            }
            runs.push((start as u32, (i - start) as u32, states[start]));
        }
        if runs.len() * std::mem::size_of::<(u32, u32, u8)>() < states.len() {
            Self::Runs(runs)
        } else {
            Self::Dense(states)
        }
    }
    pub fn at(&self, sample: usize) -> Option<State> {
        let value = match self {
            Self::None => return None,
            Self::Uniform(s) => *s,
            Self::Dense(s) => s[sample],
            Self::Runs(runs) => {
                let i = runs.partition_point(|r| r.0 as usize <= sample);
                if i == 0 || sample >= (runs[i - 1].0 as usize + runs[i - 1].1 as usize) {
                    return None;
                }
                runs[i - 1].2
            }
        };
        State::from_code(value)
    }
    pub fn validate(&self, samples: usize) -> bool {
        let valid = |value| State::from_code(value).is_some_and(|s| !s.is_passing() && s != State::RecordsAtPosition);
        match self {
            Self::None => true,
            Self::Uniform(v) => valid(*v),
            Self::Dense(v) => v.len() == samples && v.iter().all(|v| *v == PASSING || valid(*v)),
            Self::Runs(runs) => {
                let mut end = 0usize;
                runs.iter().all(|&(start, n, state)| {
                    let start = start as usize;
                    let next = start.checked_add(n as usize);
                    let good = start >= end && n != 0 && next.is_some_and(|n| n <= samples) && valid(state);
                    end = next.unwrap_or(usize::MAX);
                    good
                })
            }
        }
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct Row {
    pub missing: MissingStates,
    /// One-based record ordinals in the original VCF/BCF, excluding headers; all records count.
    pub source_records: Vec<u64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn encodings_preserve_sparse_uniform_and_dense_reasons() {
        for codes in [
            vec![PASSING; 100],
            vec![State::NoCall as u8; 100],
            (0..100)
                .map(|i| if i == 5 { State::LowQuality as u8 } else { PASSING })
                .collect(),
            (0..100)
                .map(|i| {
                    if i % 2 == 0 {
                        State::NoCall as u8
                    } else {
                        State::Filtered as u8
                    }
                })
                .collect(),
        ] {
            let states = MissingStates::from_codes(codes.clone());
            assert!(states.validate(codes.len()));
            let restored: MissingStates = serde_json::from_slice(&serde_json::to_vec(&states).unwrap()).unwrap();
            for (i, expected) in codes.iter().enumerate() {
                assert_eq!(restored.at(i), State::from_code(*expected));
            }
        }
        assert!(!MissingStates::Runs(vec![(1, 3, 6), (2, 2, 6)]).validate(10));
        assert!(!MissingStates::Uniform(200).validate(10));
        assert!(!MissingStates::Dense(vec![PASSING; 9]).validate(10));
    }
}
