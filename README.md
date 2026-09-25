# pgsum

Polygenic score calculation from single-sample gVCFs.

**Status: early development.** `compile` and `inspect` work and match the reference implementation on
11.7 million Catalog terms. `extract` and `score` are not implemented yet; the genotype rules they will use
are implemented and tested.

## What's different

Most PGS tools score a VCF and treat a site with no variant record as missing, then fill it with a mean
dosage. A gVCF says more than that: a passing reference block is a confident homozygous-reference call.
pgsum reads the gVCF directly and uses those blocks, and it never imputes. A score is reported only when
every term in the model has a usable call; otherwise it is withheld with a count of why.

## Planned usage

```sh
# Once per PGS Catalog release: compile harmonized scoring files into packs.
pgsum compile PGS000001_hmPOS_GRCh38.txt.gz --reference GRCh38.fa --out packs/

# Per sample: read genotypes at the pack sites, then score.
pgsum run --gvcf sample.g.vcf.gz --reference GRCh38.fa --pack packs/PGS000001.pgsp --out results/
```

## Scoring files and licences

pgsum does not ship any PGS weights. Download scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) yourself. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
