# pgsum

Polygenic score calculation from single-sample gVCFs.

**Status: v0.1, early development.** The full pipeline works. On GIAB HG002 it matches a reference
implementation of the same rules term for term, and sum for sum, on 11.7 million Catalog terms.

## What's different

Most PGS tools score a VCF and treat a site with no variant record as missing, then fill it with a mean
dosage. A gVCF says more than that: a passing reference block is a confident homozygous-reference call.
pgsum reads the gVCF directly and uses those blocks, and it never imputes.

Every score gets two answers:

- **strict**: a raw score only when every term in the model has a usable call (otherwise withheld, with
  the reasons);
- **partial**: the exact sum over usable terms, labelled as partial, with the share of terms and of total
  effect it covers.

Weights are summed exactly (no floating point), and every output records the digests of its inputs.

About a third of Catalog scores publish only an effect allele. They are unscorable by default; with
`--allow-inferred-other-allele`, pgsum infers the other allele from the score's own reference convention (or
the Catalog's inferred allele) and labels every result that used it. See DESIGN.md.

## Usage

```sh
# Once per PGS Catalog release: compile harmonized GRCh38 scoring files (each next to its
# <PGS_ID>.metadata.json from the Catalog REST API) into packs.
pgsum compile PGS000001_hmPOS_GRCh38.txt.gz --reference GRCh38.fa --out packs/

# Your own score: a TSV with #pgs_id= and #genome_build=GRCh38 (see DESIGN.md, "Custom scores").
pgsum compile my_score.tsv --reference GRCh38.fa --out packs/

# Per sample: read genotypes at every pack site and score, in one step ...
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/PGS000001.pgsp --out results/

# ... or in two, reusing the genotype table for more packs later.
pgsum extract --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/PGS000001.pgsp --out sample.pgsg
pgsum score --genotypes sample.pgsg --pack packs/PGS000001.pgsp --out results/ --terms

# Look inside packs and genotype tables.
pgsum inspect packs/PGS000001.pgsp [--header] [--genotypes sample.pgsg]
```

Many scores at once: `--pack` takes files or directories and can be repeated, `--pack-list` reads paths from
a file, and `--ids PGS000001,PGS000013` keeps only those scores. `extract` reads the gVCF once for all
selected packs. `score` runs packs in parallel and splits large packs into chunks across cores, so a
13-million-term score uses the whole machine. `--threads` sets the worker count for any command.

```sh
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/ --ids PGS000013,PGS000018 --out results/
```

`results/` gets `<PGS_ID>.score.json` per score, `scores.tsv` across scores, and with `--terms` a per-term
TSV.

The whole Catalog: `pgsum fetch --all --reference GRCh38.fa --out packs/` downloads and compiles all 6,990
scores (4.5 billion terms, 31.6 GB of packs); an interrupted run resumes. For many samples against the same
packs, pass `--targets-cache packs.pgst` to `extract` or `run`: the first run records which sites the packs
need, later runs skip reading every pack to find out.

On a 12-core Mac with a 36-million-record HG002 gVCF, against all 6,990 scores: extract 17 s (with the target
index; 120 s the first time) and score 95 s, so about two minutes per sample.

## Scoring files and licences

pgsum does not ship any PGS weights. Download scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) yourself. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
