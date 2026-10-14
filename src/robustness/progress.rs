//! Operational stderr events; never part of the scientific evidence or its hashes.
use std::io::Write;
use std::time::{Duration, Instant};

use clap::ValueEnum;
use serde::Serialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum ProgressFormat {
    Human,
    Json,
    Quiet,
}

#[derive(Serialize)]
struct ScoreProgress {
    pgs_id: String,
    score_index: usize,
    scores_total: usize,
    samples: usize,
    scenarios: usize,
    terms_total: usize,
    terms_processed: usize,
    baseline_terms: usize,
    excluded_terms: usize,
}

#[derive(Serialize)]
struct Event<'a> {
    schema: &'static str,
    phase: &'a str,
    elapsed_seconds: f64,
    score: Option<&'a ScoreProgress>,
    report: Option<&'a str>,
}

pub(super) struct Progress {
    format: ProgressFormat,
    interval: Duration,
    started: Instant,
    last: Duration,
    score: Option<ScoreProgress>,
}

impl Progress {
    pub(super) fn new(format: ProgressFormat, interval_seconds: u64) -> Self {
        Self {
            format,
            interval: Duration::from_secs(interval_seconds),
            started: Instant::now(),
            last: Duration::ZERO,
            score: None,
        }
    }

    pub(super) fn stage(&mut self, phase: &str) {
        self.emit(phase, None);
    }

    pub(super) fn start_score(
        &mut self,
        id: &str,
        index: usize,
        total: usize,
        samples: usize,
        scenarios: usize,
        terms: usize,
    ) {
        self.score = Some(ScoreProgress {
            pgs_id: id.into(),
            score_index: index,
            scores_total: total,
            samples,
            scenarios,
            terms_total: terms,
            terms_processed: 0,
            baseline_terms: 0,
            excluded_terms: 0,
        });
        self.stage("scoring");
    }

    /// Call only between terms. Processed includes completed exclusions, not the next term.
    pub(super) fn checkpoint(&mut self, processed: usize, eligible: usize, force: bool) {
        if self.format == ProgressFormat::Quiet {
            return;
        }
        let elapsed = self.started.elapsed();
        if !force && !due(elapsed, self.last, self.interval) {
            return;
        }
        if let Some(score) = &mut self.score {
            score.terms_processed = processed;
            score.baseline_terms = eligible;
            score.excluded_terms = processed - eligible;
        }
        self.stage("scoring");
    }

    pub(super) fn writing(&mut self) {
        self.score = None;
        self.stage("writing");
    }

    pub(super) fn complete(&mut self, report: &str) {
        self.emit("complete", Some(report));
    }

    fn emit(&mut self, phase: &str, report: Option<&str>) {
        if self.format == ProgressFormat::Quiet {
            return;
        }
        self.last = self.started.elapsed();
        let event = Event {
            schema: "pgsum-robustness-progress-v1",
            phase,
            elapsed_seconds: self.last.as_secs_f64(),
            score: self.score.as_ref(),
            report,
        };
        let line = if self.format == ProgressFormat::Json {
            serde_json::to_string(&event).expect("finite progress event")
        } else if let Some(score) = &self.score {
            format!(
                "robustness: {phase}; {} (score {}/{}; {} samples, {} scenarios); {}/{} terms processed, {} baseline, {} excluded; {:.1}s elapsed",
                score.pgs_id,
                score.score_index,
                score.scores_total,
                score.samples,
                score.scenarios,
                score.terms_processed,
                score.terms_total,
                score.baseline_terms,
                score.excluded_terms,
                event.elapsed_seconds
            )
        } else {
            format!(
                "robustness: {phase}; {:.1}s elapsed{}",
                event.elapsed_seconds,
                report.map(|path| format!("; {path}")).unwrap_or_default()
            )
        };
        // A closed log consumer must not panic or invalidate an otherwise valid experiment.
        let _ = writeln!(std::io::stderr().lock(), "{line}");
    }
}

fn due(elapsed: Duration, last: Duration, interval: Duration) -> bool {
    elapsed.saturating_sub(last) >= interval
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_is_elapsed_since_last_event_without_sleeping() {
        let seconds = Duration::from_secs;
        assert!(!due(seconds(29), seconds(0), seconds(30)));
        assert!(due(seconds(30), seconds(0), seconds(30)));
        assert!(!due(seconds(35), seconds(30), seconds(30)));
        assert!(due(seconds(60), seconds(30), seconds(30)));
        assert!(!due(seconds(5), seconds(30), seconds(30)));
    }
}
