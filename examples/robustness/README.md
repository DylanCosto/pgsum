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

The tiny sequence uses the pipeline's GRCh38 coordinate label but is entirely invented;
these are synthetic arithmetic fixtures, not human genomes or a biological PGS. Their
results establish that the packaged commands run together, not predictive accuracy,
clinical calibration, or genome-wide performance. Example input data are CC0-1.0;
the runner uses the repository's MIT OR Apache-2.0 license.

See the included [experiment guide](../../docs/robustness.md) for real input preparation,
scope restrictions, and interpretation.
