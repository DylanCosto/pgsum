//! Self-contained, script-free HTML summary; detailed evidence stays in linked machine-readable files.
use super::Report;
use crate::{Error, Result};
use std::fmt::Write as _;
use std::path::Path;

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
fn number(v: Option<f64>) -> String {
    v.map(|x| {
        if x != 0.0 && (x.abs() < 1e-4 || x.abs() >= 1e6) {
            format!("{x:.5e}")
        } else {
            format!("{x:.5}")
        }
    })
    .unwrap_or_else(|| "Unavailable".into())
}

pub(super) fn write(path: &Path, report: &Report) -> Result<()> {
    let mut html = String::from(
        r##"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PGS robustness experiment</title><style>
:root{color-scheme:light;--ink:#172e38;--muted:#526773;--line:#d5e2e7;--accent:#006b66}*{box-sizing:border-box}body{font:16px/1.55 system-ui,sans-serif;color:var(--ink);background:#f3f7f8;margin:0}main{max-width:1280px;margin:0 auto;padding:40px 28px}header{border-bottom:3px solid var(--accent);padding-bottom:20px;margin-bottom:28px}.eyebrow{color:var(--accent);font-weight:700;letter-spacing:.12em;font-size:12px;text-transform:uppercase}h1{font-size:36px;line-height:1.2;margin:8px 0}h2{margin:0 0 8px}p{max-width:86ch}.muted,small{color:var(--muted)}section{background:white;border:1px solid var(--line);border-radius:12px;padding:24px;margin:22px 0}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:14px}th,td{text-align:right;padding:12px;border-bottom:1px solid var(--line);vertical-align:top}th:first-child,td:first-child{text-align:left}th{background:#edf5f5;color:var(--accent);white-space:nowrap}a{color:#005b92}nav{display:flex;gap:20px;flex-wrap:wrap}.pill{display:inline-block;background:#edf5f5;border-radius:5px;padding:4px 9px;font-size:13px;margin:2px}details{margin-top:18px}summary{cursor:pointer;font-weight:650}li{margin:8px 0}code{font-size:13px;overflow-wrap:anywhere}@media(max-width:600px){main{padding:22px 14px}section{padding:16px}h1{font-size:28px}}
</style></head><body><main><header><div class="eyebrow">pgsum / Research diagnostics</div><h1>How stable are these polygenic scores?</h1><p>Controlled changes to genotype availability and analysis inputs, measured against a frozen available-data baseline. These results describe score reconstruction and numerical rank stability; they do not measure disease prediction.</p><nav><a href="samples.tsv">Sample results</a><a href="top-changes.jsonl">Variant explanations</a><a href="references.tsv">Reference comparisons</a><a href="regions.tsv">Region contributions</a><a href="report.json">Full report &amp; provenance</a></nav></header>"##,
    );
    html.push_str("<section><h2>Read the experiment correctly</h2><ul>");
    for note in &report.notes {
        let _ = write!(html, "<li>{}</li>", escape(note));
    }
    html.push_str("</ul><p><strong>Score errors are in each model’s own units.</strong> Compare scenarios within a model. Percentile shifts are percentage points against the unchanged baseline distribution; Spearman correlation measures ordering after both cohorts are ranked separately.</p></section>");
    for score in &report.scores {
        let _ = write!(
            html,
            "<section><h2>{}</h2><div class=\"pill\">{} / {} published-file terms in baseline ({:.2}%)</div><div class=\"pill\">Weight coverage: {}</div><p class=\"muted\">The same eligible terms are used for every scenario and selected reference set. Scope: {}.</p>",
            escape(&score.pgs_id),
            score.baseline_terms,
            score.total_model_terms,
            100.0 * score.baseline_term_fraction,
            score
                .baseline_weight_fraction
                .map(|x| format!("{:.2}%", 100.0 * x))
                .unwrap_or_else(|| "unavailable".into()),
            if score.baseline_terms == score.total_model_terms
                && score.inventory_consistent
                && score.matches_publication != Some(false)
            {
                "all terms in this compiled model"
            } else {
                "a conditional subset of this model"
            }
        );
        html.push_str("<div class=\"scroll\"><table><thead><tr><th>Scenario</th><th>Samples</th><th>Mean absolute error</th><th>RMSE</th><th>Spearman</th><th>Mean absolute percentile shift</th><th>Largest percentile shift</th><th>Masked / filled calls</th></tr></thead><tbody>");
        for s in &score.scenarios {
            let m = &s.metrics;
            let _ = write!(
                html,
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{} / {}</td></tr>",
                escape(&s.id),
                m.n,
                number(m.mae),
                number(m.rmse),
                number(m.spearman),
                number(m.mean_abs_percentile_shift),
                number(m.max_abs_percentile_shift),
                s.masked_observations,
                s.filled_observations
            );
        }
        html.push_str("</tbody></table></div><details><summary>Cutoff crossings and exclusions</summary><p>Cutoff membership means percentile ≥ cutoff. Counts are paired changes against the fixed baseline distribution.</p><div class=\"scroll\"><table><thead><tr><th>Scenario</th><th>Cutoff</th><th>Up</th><th>Down</th><th>Unchanged</th><th>Evaluable</th></tr></thead><tbody>");
        for s in &score.scenarios {
            for c in &s.metrics.crossings {
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                    escape(&s.id),
                    c.cutoff,
                    c.up,
                    c.down,
                    c.unchanged,
                    c.evaluable
                );
            }
        }
        html.push_str("</tbody></table></div><ul>");
        for (reason, n) in &score.excluded_terms {
            let _ = write!(html, "<li>{}: {}</li>", escape(reason), n);
        }
        html.push_str("</ul></details>");
        if score.scenarios.iter().any(|s| !s.reference_comparisons.is_empty()) {
            html.push_str("<details><summary>Frozen reference comparisons</summary><p>Each row uses the same baseline terms. Mean percentiles summarize numerical placement, not population risk. Shared sample identities indicate that the reference is not independent of the target cohort.</p><div class=\"scroll\"><table><thead><tr><th>Scenario / reference</th><th>Reference samples</th><th>Shared identities</th><th>Mean baseline percentile</th><th>Mean scenario percentile</th><th>Mean absolute shift</th></tr></thead><tbody>");
            for s in &score.scenarios {
                for r in &s.reference_comparisons {
                    let _ = write!(
                        html,
                        "<tr><td>{} / {}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                        escape(&s.id),
                        escape(&r.id),
                        r.reference_samples,
                        r.shared_sample_ids,
                        number(r.mean_baseline_percentile),
                        number(r.mean_scenario_percentile),
                        number(r.mean_absolute_shift)
                    );
                }
            }
            html.push_str("</tbody></table></div></details>");
        }
        html.push_str("</section>");
    }
    html.push_str("<section><h2>Experiment choices</h2><ul>");
    for s in &report.experiment.scenarios {
        let settings = serde_json::to_string(s).map_err(|e| Error::Invalid(e.to_string()))?;
        let _ = write!(
            html,
            "<li><strong>{}</strong><br><code>{}</code></li>",
            escape(&s.id),
            escape(&settings)
        );
    }
    html.push_str("</ul><p>BED intervals use zero-based, half-open coordinates and test the variant anchor position. Random masks are deterministic for the declared seed, position, and (for sample-site masks) sample identity. Quality or dosage-policy comparisons use separately extracted genotype tables.</p></section><footer class=\"muted\">Generated by pgSum. An experiment is complete only when report.json exists and its recorded output hashes match these files.</footer></main></body></html>");
    std::fs::write(path, html).map_err(Error::io(path))
}

#[cfg(test)]
mod tests {
    #[test]
    fn escaped_html() {
        assert_eq!(super::escape("<&\"'>"), "&lt;&amp;&quot;&#39;&gt;");
    }
    #[test]
    fn small_errors_remain_visible() {
        assert_eq!(super::number(Some(1e-8)), "1.00000e-8");
    }
}
