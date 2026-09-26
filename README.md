# pgsum

Polygenic score calculation from gVCFs, exact and without imputation.

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

By default pgsum follows a conservative reference implementation exactly. Opt-ins widen what can be scored,
each validated on real data and labelled in every result that uses it (see DESIGN.md):

| Option | Scores terms that… |
|---|---|
| `--allow-inferred-other-allele` | publish only an effect allele (about 30% of Catalog terms) |
| `--allow-inferred-palindromes` | are A/T or C/G SNVs, in scores whose other SNVs are on the forward strand |
| `--allow-inferred-indels` | are indels or multi-base variants (needs `compile --public-variants`) |
| `--accept-informational-descriptions` | carry a `variant_description` that is only an annotation |
| `extract --haploid-xy-as-homozygous` | sit on haploid chrX/chrY calls from callers that write them |

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
TSV. With `--bundle`, every result goes into one `results.jsonl.zst` (one JSON object per line) instead of a
file per score: for the whole Catalog that is 4 MB rather than 6,990 small files, which on an exFAT drive with
1 MB allocation blocks take 14 GB.

### Inputs

pgsum reads gVCFs from DeepVariant (long- and short-read), DRAGEN and GATK HaplotypeCaller, and plain VCFs,
on GRCh38:

- bgzipped (`.vcf.gz`), plain gzip or uncompressed; BCF is not read (convert with `bcftools view -Oz`);
- contigs named `chr1` or `1` (and `chrM` or `MT`), in the VCF and in the reference FASTA alike;
- one sample, or several with `--sample NAME`;
- a `.tbi` or `.csi` index next to a bgzipped file lets `extract` read only the parts the scores need, so a
  run with a few scores takes a fraction of a second (`tabix -p vcf sample.g.vcf.gz` creates one). pgsum
  decides per run (`--scan auto`, the default) and the table is identical either way.

A VCF whose `##contig` lengths differ from the reference (for example GRCh37) is refused. A plain VCF has no
reference blocks, so sites without a record are `unknown_no_record` rather than assumed reference. DRAGEN
writes male chrX and chrY as haploid; pass `--haploid-xy-as-homozygous` to `extract` or `run`.

pgsum records the SHA-256 of every input. Hashing a 3 GB reference takes a second, so digests are remembered
in `~/.cache/pgsum/digests.tsv` (or `$XDG_CACHE_HOME/pgsum`, or `$PGSUM_CACHE_DIR`), keyed by path, size,
modification time and inode; `PGSUM_NO_DIGEST_CACHE=1` turns this off.

The whole Catalog: `pgsum fetch --all --reference GRCh38.fa --out packs/` downloads and compiles all 6,990
scores (4.5 billion terms, 31.6 GB of packs); an interrupted run resumes. For many samples against the same
packs, pass `--targets-cache packs.pgst` to `extract` or `run`: the first run records which sites the packs
need, later runs skip reading every pack to find out.

### Speed

On a 12-core Mac with the packs on a USB SSD and a 36-million-record HG002 gVCF:

| Run | Time |
|---|---|
| `run`, one score (77 terms), indexed gVCF | 0.4 s |
| `run`, four small scores (588 targets), indexed gVCF | 0.8 s |
| `extract`, all 6,990 scores (43.3M targets, with the target index) | 12 s |
| `score`, all 6,990 scores (4.5 billion terms) | 141 s |
| `score`, three scores against a saved whole-Catalog genotype table | 0.8 s |

Without an index, `extract` reads the whole file on all cores: 5 s for a 154-million-record GATK gVCF.

## Scoring files and licences

pgsum does not ship any PGS weights. Download scoring files from the
[PGS Catalog](https://www.pgscatalog.org/) yourself. Individual scores carry their own licence terms, and
pgsum keeps each score's licence field in its output.

## Design

See [DESIGN.md](DESIGN.md) for the pipeline, genotype rules, completeness rules and output format.

## Licence

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at
your option.
