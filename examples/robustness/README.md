# Offline robustness example

This miniature example includes every input and uses no network connection. It needs
the `pgsum` executable for your platform and Python 3.8 or newer (standard library only).
From an unpacked release directory:

```sh
python3 examples/robustness/run.py --binary ./pgsum --out example-output
```

Use a new output directory. The runner works from any current directory when its paths
are supplied explicitly, and leaves the example inputs unchanged. In a source checkout,
build the executable first and use `--binary target/release/pgsum`.

The runner compiles the three-term model, extracts four synthetic targets and four
separate synthetic reference samples, builds the frozen frequency table, and runs four
scenarios: unchanged control, omission of the last site, omission with mean filling,
and all sites missing. It then checks exact scores/differences, counts, percentiles,
reference comparisons, regional attribution, and all reported output hashes against
the hand-calculated expectations in `expected.json`. Successful verification writes
`example-output/verification.json`; the report is `example-output/results/report.html`.

The baseline scores for A-D are `3, -1.5, 3, 1.5`. Dropping the last site gives
`-1, -1.5, 1, -0.5`. Its ALT frequency among R1-R4 is `6/8 = 0.75`, so filling adds
`2 × 2 × 0.75 = 3` per sample and gives `2, 1.5, 4, 2.5`. This deliberately includes
a masked zero contribution and tied baseline scores. Samples A-D never contribute to
reference frequencies.

## Paired missingness study

To try the complete planning and comparison workflow, use Python 3.9 or newer:

```sh
python3 examples/robustness/run.py --binary ./pgsum --study --out study-example-output
```

This also plans omission versus frequency filling at 0%, 10%, 20% and 100% missingness,
with three reproducible sample-site masking seeds. The 26 scenarios include two controls
and produce 104 participant/scenario rows. It runs `robustness-plan`, `robustness` and
`robustness-summarize`, then checks score metrics, cutoff transitions and frozen-reference
shifts independently from the participant output files. The 100% omission boundary checks
that undefined correlations stay unavailable.

Open `study-example-output/study-summary/summary.html` to compare policies and seed
variation. JSON and TSV files retain full precision; `study-verification.json` records
the independent checks. `policies.json`, `study-plan/` and `commands.json` show the complete
reproducible workflow. The original hand-calculated example still runs and is checked first.

The tiny sequence uses the pipeline's GRCh38 coordinate label but is entirely invented;
these are synthetic arithmetic fixtures, not human genomes or a biological PGS. Their
results establish that the packaged commands run together, not predictive accuracy,
clinical calibration, or genome-wide performance. Example input data are CC0-1.0;
the runner uses the repository's MIT OR Apache-2.0 license.

See the included [experiment guide](../../docs/robustness.md) for real input preparation,
scope restrictions, and interpretation.

The paired study also runs `robustness-groups` with one synthetic singleton and one three-person
group. Open `study-groups/groups.html` to inspect score errors, whole-cohort cutoff transitions
and frozen-reference shifts by group. `group-verification.json` records the independent check;
these artificial labels have no demographic or biological meaning.

The same mapping is passed to `robustness-summarize --groups` to compare policies by rate and
seed within each group. Open `group-study/summary.html`; `group-study-verification.json`
records independent checks against the native participant rows. Singleton correlations remain
unavailable while percentile shifts retain the full four-person target baseline.
