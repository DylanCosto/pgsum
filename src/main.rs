use clap::{Args, Parser, Subcommand};
use pgsum::compile::{compile_file, reference_identity};
use pgsum::extract::{Options as ExtractOptions, ScanMode};
use pgsum::genotypes::GenotypeTable;
use pgsum::pack::Pack;
use pgsum::reference::Reference;
use pgsum::{Error, Result};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Polygenic score calculation from single-sample gVCFs.
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    /// Worker threads for every step (default: all cores).
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Args)]
struct FetchArgs {
    /// These PGS IDs (comma-separated or repeated).
    #[arg(long, value_delimiter = ',', required_unless_present = "all", conflicts_with = "all")]
    ids: Vec<String>,
    /// Every score in the Catalog.
    #[arg(long)]
    all: bool,
    /// GRCh38 reference FASTA (with `.fai`).
    #[arg(long)]
    reference: PathBuf,
    /// Directory to write packs into.
    #[arg(long)]
    out: PathBuf,
    /// Where downloads go while they are compiled (default: `<out>/downloads`).
    #[arg(long)]
    downloads: Option<PathBuf>,
    /// Keep scoring files and Catalog records after compiling.
    #[arg(long)]
    keep_downloads: bool,
    /// Skip scores with more variants than this.
    #[arg(long)]
    max_variants: Option<u64>,
    /// Scores downloaded and compiled at once.
    #[arg(long, default_value_t = 4)]
    jobs: usize,
    /// Public variant set (CHROM POS REF ALT) used to orient indels; see DESIGN.md.
    #[arg(long)]
    public_variants: Option<PathBuf>,
}

/// Which packs to use. `--pack` takes files or directories (every `*.pgsp` inside) and can be repeated;
/// `--pack-list` adds paths from a file; `--ids` keeps only the listed scores.
#[derive(Args, Clone)]
struct PackSelection {
    /// A pack, or a directory of packs. Repeatable.
    #[arg(long = "pack")]
    packs: Vec<PathBuf>,
    /// A file listing pack paths or directories, one per line (`#` comments allowed).
    #[arg(long)]
    pack_list: Option<PathBuf>,
    /// Keep only these score IDs (comma-separated or repeated), e.g. `PGS000001,PGS000013,MY_SCORE`.
    #[arg(long, value_delimiter = ',')]
    ids: Vec<String>,
}

impl PackSelection {
    fn resolve(&self) -> Result<Vec<PathBuf>> {
        let mut roots = self.packs.clone();
        if let Some(list) = &self.pack_list {
            let text = std::fs::read_to_string(list).map_err(Error::io(list))?;
            let base = list.parent().unwrap_or(Path::new("."));
            for line in text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
            {
                roots.push(base.join(line));
            }
        }
        let mut paths = Vec::new();
        for root in roots {
            if root.is_dir() {
                for entry in std::fs::read_dir(&root).map_err(Error::io(&root))? {
                    let path = entry.map_err(Error::io(&root))?.path();
                    if path.extension().is_some_and(|e| e == pgsum::pack::EXTENSION) && !is_hidden(&path) {
                        paths.push(path);
                    }
                }
            } else {
                paths.push(root);
            }
        }
        if !self.ids.is_empty() {
            for id in &self.ids {
                if !pgsum::scoring_file::is_pgs_id(id) && !pgsum::scoring_file::is_custom_id(id) {
                    return Err(Error::Invalid(format!(
                        "{id:?} is not a PGS Catalog or custom score ID"
                    )));
                }
            }
            paths.retain(|p| {
                p.file_stem()
                    .and_then(|s| s.to_str())
                    .is_some_and(|s| self.ids.iter().any(|id| id == s))
            });
            let found: Vec<_> = paths.iter().filter_map(|p| p.file_stem()?.to_str()).collect();
            let missing: Vec<_> = self.ids.iter().filter(|id| !found.contains(&id.as_str())).collect();
            if !missing.is_empty() {
                return Err(Error::Invalid(format!("no pack for {missing:?}")));
            }
        }
        paths.sort();
        paths.dedup();
        if paths.is_empty() {
            return Err(Error::Invalid("no packs selected; use --pack or --pack-list".into()));
        }
        Ok(paths)
    }
}

#[derive(Subcommand)]
enum Command {
    /// Compile PGS Catalog harmonized scoring files into packs.
    Compile {
        /// Scoring files: PGS Catalog harmonized files (`*_hmPOS_GRCh38.txt.gz`, or directories of them, each with
        /// its `<PGS_ID>.metadata.json` alongside), or custom scoring files (see DESIGN.md, "Custom scores").
        #[arg(required = true)]
        scoring_files: Vec<PathBuf>,
        /// GRCh38 reference FASTA (with `.fai`), used to resolve SNV orientation.
        #[arg(long)]
        reference: PathBuf,
        /// Directory to write packs into.
        #[arg(long)]
        out: PathBuf,
        /// Public variant set (CHROM POS REF ALT) used to orient indels; see DESIGN.md.
        #[arg(long)]
        public_variants: Option<PathBuf>,
    },
    /// Download PGS Catalog scores and compile each into a pack, deleting downloads as it goes.
    ///
    /// Scores whose pack already exists with an identical Catalog record are skipped, so an interrupted
    /// run can be restarted. A log of every score is written to `<out>/fetch.tsv`.
    Fetch(FetchArgs),
    /// Print a pack's terms as TSV, or its header as JSON.
    Inspect {
        /// A pack (prints its terms) or a genotype table (prints its targets and calls).
        path: PathBuf,
        /// Print the header instead of the terms or targets.
        #[arg(long)]
        header: bool,
        /// Add each term's status, call state and effect dosage from this genotype table.
        #[arg(long)]
        genotypes: Option<PathBuf>,
        /// With `--genotypes`: allow inferred other alleles, as `score` does.
        #[arg(long)]
        allow_inferred_other_allele: bool,
        /// With `--genotypes`: accept informational variant descriptions, as `score` does.
        #[arg(long)]
        accept_informational_descriptions: bool,
        /// With `--genotypes`: allow inferred palindromes, as `score` does.
        #[arg(long)]
        allow_inferred_palindromes: bool,
        /// With `--genotypes`: allow inferred indels, as `score` does.
        #[arg(long)]
        allow_inferred_indels: bool,
    },
    /// Print a pack's per-term evidence rows (JSON lines): status, call, effect dosage, contribution and the gVCF
    /// records at the term's position. The table must be extracted with `--term-positions`.
    Evidence {
        pack: PathBuf,
        #[arg(long)]
        genotypes: PathBuf,
        /// GRCh38 reference FASTA (with `.fai`).
        #[arg(long)]
        reference: PathBuf,
        /// Policy ID written into each call (default: the genotype table's policy ID).
        #[arg(long)]
        call_policy: Option<String>,
        /// Start after this many terms (to resume).
        #[arg(long, default_value_t = 0)]
        after: usize,
    },
    /// Read genotypes from a gVCF at every site the packs need.
    Extract {
        /// Bgzipped single-sample gVCF with a tabix or CSI index.
        #[arg(long)]
        gvcf: PathBuf,
        /// GRCh38 reference FASTA (with `.fai`).
        #[arg(long)]
        reference: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        /// Output genotype table (`.pgsg`).
        #[arg(long)]
        out: PathBuf,
        /// Target index to reuse when it matches the selected packs, or to create (`.pgst`).
        #[arg(long)]
        targets_cache: Option<PathBuf>,
        /// Read haploid chrX/chrY calls (`1`) as homozygous (`1/1`), as DeepVariant writes male chrX.
        #[arg(long)]
        haploid_xy_as_homozygous: bool,
        /// How to read the gVCF: `auto` uses its .tbi/.csi index when the targets need under half of the
        /// file, `full` always reads it whole, `indexed` requires the index.
        #[arg(long, value_enum, default_value_t = ScanMode::Auto)]
        scan: ScanMode,
        /// The sample to read from a multi-sample VCF.
        #[arg(long)]
        sample: Option<String>,
        /// Read every sample of a multi-sample VCF into a cohort file (`.pgsc`), in one pass; `score` then
        /// scores every sample.
        #[arg(long, conflicts_with = "sample")]
        all_samples: bool,
        /// Accept calls that report neither depth nor GQ, on GT and FILTER alone (genotype-only VCFs: imputed,
        /// array or joint-called data). Recorded in the table's policy ID.
        #[arg(long)]
        accept_missing_quality: bool,
        /// Leave out records whose ALTs are all structural-variant symbols (`<DEL>`, `<INV>`, …, breakends), which
        /// otherwise make every target they span ambiguous; for panels that carry SVs, such as the 30× 1000
        /// Genomes release. Recorded in the policy ID.
        #[arg(long)]
        skip_structural_alleles: bool,
        /// Read records that start at the same position (a multi-allelic site split into one record per ALT) as
        /// one multi-allelic record, instead of calling the site ambiguous; for panels written that way, such as
        /// the 30× 1000 Genomes release. Recorded in the policy ID.
        #[arg(long)]
        merge_split_records: bool,
        /// Also keep the gVCF records at the position of every term with one, including terms that cannot be
        /// scored (see `evidence`).
        #[arg(long)]
        term_positions: bool,
    },
    /// Score packs against an extracted genotype table.
    Score {
        /// Genotype table from `extract`, or a cohort file from `extract --all-samples`.
        #[arg(long)]
        genotypes: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        /// Output directory for one JSON result per score.
        #[arg(long)]
        out: PathBuf,
        /// Also write a per-term TSV for each score.
        #[arg(long)]
        terms: bool,
        /// Write every result to one `results.jsonl.zst` (one JSON object per line, sorted by score ID) instead
        /// of one file per score. Much smaller on drives with large allocation blocks (e.g. exFAT).
        #[arg(long)]
        bundle: bool,
        /// Also score terms without an author other allele, using the orientation their pack inferred.
        #[arg(long)]
        allow_inferred_other_allele: bool,
        /// Also score terms whose variant_description is informational (fine-mapping statistics, variant IDs).
        #[arg(long)]
        accept_informational_descriptions: bool,
        /// Also score palindromic SNVs on the forward strand when the score's other SNVs are all there.
        #[arg(long)]
        allow_inferred_palindromes: bool,
        /// Also score indels and multi-base terms oriented by reference fit or a public variant set.
        #[arg(long)]
        allow_inferred_indels: bool,
        #[command(flatten)]
        reference: ReferenceArgs,
    },
    /// Extract then score in one step.
    Run {
        #[arg(long)]
        gvcf: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[command(flatten)]
        packs: PackSelection,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        terms: bool,
        /// Write every result to one `results.jsonl.zst` (one JSON object per line, sorted by score ID) instead
        /// of one file per score. Much smaller on drives with large allocation blocks (e.g. exFAT).
        #[arg(long)]
        bundle: bool,
        /// Target index to reuse when it matches the selected packs, or to create (`.pgst`).
        #[arg(long)]
        targets_cache: Option<PathBuf>,
        /// Read haploid chrX/chrY calls (`1`) as homozygous (`1/1`), as DeepVariant writes male chrX.
        #[arg(long)]
        haploid_xy_as_homozygous: bool,
        /// How to read the gVCF: `auto` uses its .tbi/.csi index when the targets need under half of the
        /// file, `full` always reads it whole, `indexed` requires the index.
        #[arg(long, value_enum, default_value_t = ScanMode::Auto)]
        scan: ScanMode,
        /// The sample to read from a multi-sample VCF.
        #[arg(long)]
        sample: Option<String>,
        /// Accept calls that report neither depth nor GQ, on GT and FILTER alone (genotype-only VCFs: imputed,
        /// array or joint-called data). Recorded in the table's policy ID.
        #[arg(long)]
        accept_missing_quality: bool,
        /// Leave out records whose ALTs are all structural-variant symbols (`<DEL>`, `<INV>`, …, breakends), which
        /// otherwise make every target they span ambiguous; for panels that carry SVs, such as the 30× 1000
        /// Genomes release. Recorded in the policy ID.
        #[arg(long)]
        skip_structural_alleles: bool,
        /// Read records that start at the same position (a multi-allelic site split into one record per ALT) as
        /// one multi-allelic record, instead of calling the site ambiguous; for panels written that way, such as
        /// the 30× 1000 Genomes release. Recorded in the policy ID.
        #[arg(long)]
        merge_split_records: bool,
        /// Also score terms without an author other allele, using the orientation their pack inferred.
        #[arg(long)]
        allow_inferred_other_allele: bool,
        /// Also score terms whose variant_description is informational (fine-mapping statistics, variant IDs).
        #[arg(long)]
        accept_informational_descriptions: bool,
        /// Also score palindromic SNVs on the forward strand when the score's other SNVs are all there.
        #[arg(long)]
        allow_inferred_palindromes: bool,
        /// Also score indels and multi-base terms oriented by reference fit or a public variant set.
        #[arg(long)]
        allow_inferred_indels: bool,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let threads = threads(cli.threads);
    if let Err(e) = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global() {
        eprintln!("pgsum: {e}");
        return ExitCode::FAILURE;
    }
    let result: Result<()> = match cli.command {
        Command::Compile {
            scoring_files,
            reference,
            out,
            public_variants,
        } => compile(&scoring_files, &reference, &out, public_variants.as_deref()),
        Command::Inspect {
            path,
            header,
            genotypes,
            allow_inferred_other_allele,
            accept_informational_descriptions,
            allow_inferred_palindromes,
            allow_inferred_indels,
        } => inspect(
            &path,
            header,
            genotypes.as_deref(),
            &pgsum::score::Options {
                allow_inferred_other_allele,
                accept_informational_descriptions,
                allow_inferred_palindromes,
                allow_inferred_indels,
            },
        ),
        Command::Evidence {
            pack,
            genotypes,
            reference,
            call_policy,
            after,
        } => evidence(&pack, &genotypes, &reference, call_policy.as_deref(), after),
        Command::Fetch(args) => fetch(&args),
        Command::Extract {
            gvcf,
            reference,
            packs,
            out,
            targets_cache,
            haploid_xy_as_homozygous,
            scan,
            sample,
            accept_missing_quality,
            skip_structural_alleles,
            merge_split_records,
            term_positions,
            all_samples,
        } => packs.resolve().and_then(|packs| {
            if all_samples {
                let options = pgsum::cohort::CohortOptions {
                    targets_cache: targets_cache.as_deref(),
                    haploid_xy_as_homozygous,
                    accept_missing_quality,
                    skip_structural_alleles,
                    merge_split_records,
                    threads,
                };
                return extract_cohort(&gvcf, &reference, &packs, &out, &options);
            }
            let options = ExtractOptions {
                targets_cache: targets_cache.as_deref(),
                haploid_xy_as_homozygous,
                accept_missing_quality,
                skip_structural_alleles,
                merge_split_records,
                term_positions,
                sample: sample.as_deref(),
                scan,
                threads,
            };
            extract(&gvcf, &reference, &packs, &out, &options)
        }),
        Command::Score {
            genotypes,
            packs,
            out,
            terms,
            bundle,
            allow_inferred_other_allele,
            accept_informational_descriptions,
            allow_inferred_palindromes,
            allow_inferred_indels,
            reference,
        } => packs.resolve().and_then(|packs| {
            let options = pgsum::score::Options {
                allow_inferred_other_allele,
                accept_informational_descriptions,
                allow_inferred_palindromes,
                allow_inferred_indels,
            };
            if is_cohort(&genotypes)? {
                return score_cohort(&genotypes, &packs, &out, &options);
            }
            let table = GenotypeTable::open(&genotypes)?;
            let panel = open_panel(&reference, &table)?;
            score(&table, &packs, &out, terms, bundle, &options, panel.as_ref())
        }),
        Command::Run {
            gvcf,
            reference,
            packs,
            out,
            terms,
            bundle,
            targets_cache,
            haploid_xy_as_homozygous,
            scan,
            sample,
            accept_missing_quality,
            skip_structural_alleles,
            merge_split_records,
            allow_inferred_other_allele,
            accept_informational_descriptions,
            allow_inferred_palindromes,
            allow_inferred_indels,
        } => packs.resolve().and_then(|packs| {
            let options = pgsum::score::Options {
                allow_inferred_other_allele,
                accept_informational_descriptions,
                allow_inferred_palindromes,
                allow_inferred_indels,
            };
            std::fs::create_dir_all(&out).map_err(Error::io(&out))?;
            let table_path = out.join("genotypes.pgsg");
            let extract_options = ExtractOptions {
                targets_cache: targets_cache.as_deref(),
                haploid_xy_as_homozygous,
                accept_missing_quality,
                skip_structural_alleles,
                merge_split_records,
                term_positions: false,
                sample: sample.as_deref(),
                scan,
                threads,
            };
            extract(&gvcf, &reference, &packs, &table_path, &extract_options)?;
            score(
                &GenotypeTable::open(&table_path)?,
                &packs,
                &out,
                terms,
                bundle,
                &options,
                None,
            )
        }),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        // The reader of stdout went away (e.g. `pgsum inspect … | head`): not an error.
        Err(Error::Io { source, .. }) if source.kind() == std::io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("pgsum: {e}");
            ExitCode::FAILURE
        }
    }
}

/// A reference panel to place scores in (`score --reference`).
#[derive(clap::Args, Clone, Debug, Default)]
struct ReferenceArgs {
    /// Place each score among this panel's scores over the same terms: a cohort file from
    /// `extract --all-samples` with the same packs (for example 1000 Genomes).
    #[arg(long = "reference-panel", requires = "reference_groups")]
    reference_panel: Option<PathBuf>,
    /// TSV naming each panel sample's group: a header with `sample` and the group column.
    #[arg(long = "reference-groups")]
    reference_groups: Option<PathBuf>,
    /// The group column of `--reference-groups`.
    #[arg(long = "reference-group-column", default_value = "super_pop")]
    reference_group_column: String,
}

fn threads(requested: Option<usize>) -> usize {
    requested
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
        .max(1)
}

/// Dot files, including the `._name` metadata files macOS writes next to files on exFAT and FAT drives.
fn is_hidden(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.starts_with('.'))
}

/// Scoring files given directly, plus every `*_hmPOS_GRCh38.txt.gz` in given directories.
fn scoring_file_paths(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for input in inputs {
        if input.is_dir() {
            for entry in std::fs::read_dir(input).map_err(Error::io(input))? {
                let path = entry.map_err(Error::io(input))?.path();
                if path.to_str().is_some_and(|p| p.ends_with("_hmPOS_GRCh38.txt.gz")) && !is_hidden(&path) {
                    paths.push(path);
                }
            }
        } else {
            paths.push(input.clone());
        }
    }
    paths.sort();
    paths.dedup();
    Ok(paths)
}

/// Load a public variant set, reporting its size.
fn public_variants(path: Option<&Path>) -> Result<Option<pgsum::public::PublicVariants>> {
    let Some(path) = path else { return Ok(None) };
    let p = pgsum::public::PublicVariants::open(path)?;
    eprintln!("{}: {} public variant records", path.display(), p.len());
    Ok(Some(p))
}

fn compile(inputs: &[PathBuf], reference: &Path, out: &Path, public: Option<&Path>) -> Result<()> {
    let public = public_variants(public)?;
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let scoring_files = scoring_file_paths(inputs)?;
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let failures: Vec<String> = scoring_files
        .par_iter()
        .filter_map(|path| {
            let started = std::time::Instant::now();
            let result = compile_file(path, None, &reference, &identity, public.as_ref(), out);
            match result {
                Ok(h) => {
                    eprintln!(
                        "{}: {} terms, {} resolved, inventory {} ({:.1}s)",
                        h.pgs_id,
                        h.inventory.actual_terms,
                        h.counts.orientation.get("resolved").copied().unwrap_or(0),
                        if h.inventory.consistent {
                            "consistent"
                        } else {
                            "INCONSISTENT"
                        },
                        started.elapsed().as_secs_f64()
                    );
                    None
                }
                Err(e) => Some(e.to_string()),
            }
        })
        .collect();
    if failures.is_empty() {
        Ok(())
    } else {
        Err(Error::Invalid(failures.join("\n")))
    }
}

fn fetch(args: &FetchArgs) -> Result<()> {
    let public = public_variants(args.public_variants.as_deref())?;
    use pgsum::fetch::Outcome;
    let started = std::time::Instant::now();
    let reference = Reference::open(&args.reference)?;
    let identity = reference_identity(&reference)?;
    let agent = pgsum::fetch::agent();
    let records = pgsum::fetch::catalog_records(&agent, &args.ids)?;
    let variants: u64 = records.iter().filter_map(|r| r["variants_number"].as_u64()).sum();
    eprintln!("{} scores, {variants} variants in the Catalog records", records.len());
    let out = args.out.as_path();
    let downloads = args.downloads.clone().unwrap_or_else(|| out.join("downloads"));
    let options = pgsum::fetch::Options {
        reference: &reference,
        identity: &identity,
        out,
        downloads: &downloads,
        keep_downloads: args.keep_downloads,
        max_variants: args.max_variants,
        public: public.as_ref(),
    };
    let done = std::sync::atomic::AtomicUsize::new(0);
    let total = records.len();
    let report = |l: &pgsum::fetch::Logged| {
        let n = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        let what = match &l.outcome {
            Outcome::Compiled { terms, .. } => format!("compiled, {terms} terms"),
            Outcome::Unchanged => "unchanged".into(),
            Outcome::Skipped(why) => format!("skipped: {why}"),
            Outcome::Failed(why) => format!("FAILED: {why}"),
        };
        eprintln!("[{n}/{total}] {}: {what} ({:.1}s)", l.id, l.seconds);
    };
    let mut log = pgsum::fetch::fetch_all(&agent, &records, &options, args.jobs, &report)?;
    log.sort_by(|a, b| a.id.cmp(&b.id));
    let mut tsv = String::from("pgs_id\toutcome\tterms\tsource_bytes\tseconds\tdetail\n");
    let mut failed = 0;
    for l in &log {
        let (outcome, terms, bytes, detail) = match &l.outcome {
            Outcome::Compiled { terms, source_bytes } => ("compiled", terms.to_string(), source_bytes.to_string(), ""),
            Outcome::Unchanged => ("unchanged", String::new(), String::new(), ""),
            Outcome::Skipped(why) => ("skipped", String::new(), String::new(), why.as_str()),
            Outcome::Failed(why) => {
                failed += 1;
                ("failed", String::new(), String::new(), why.as_str())
            }
        };
        tsv.push_str(&format!(
            "{}\t{outcome}\t{terms}\t{bytes}\t{:.1}\t{}\n",
            l.id,
            l.seconds,
            detail.replace(['\t', '\n'], " ")
        ));
    }
    let log_path = out.join("fetch.tsv");
    std::fs::write(&log_path, tsv).map_err(Error::io(&log_path))?;
    eprintln!(
        "done in {:.0}s; log in {}",
        started.elapsed().as_secs_f64(),
        log_path.display()
    );
    if failed > 0 {
        return Err(Error::Invalid(format!(
            "{failed} of {total} scores failed; see {}",
            log_path.display()
        )));
    }
    Ok(())
}

fn extract(gvcf: &Path, reference: &Path, packs: &[PathBuf], out: &Path, options: &ExtractOptions) -> Result<()> {
    let started = std::time::Instant::now();
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let loaded = started.elapsed().as_secs_f64();
    let (table, timings) = pgsum::extract::extract(gvcf, &reference, &identity, packs, options)?;
    let t = std::time::Instant::now();
    table.write(out)?;
    eprintln!(
        "  timings: reference {loaded:.1}s, targets {:.1}s{}, scan {:.1}s, assess {:.1}s, table {:.1}s, write {:.1}s",
        timings.targets_s,
        if timings.targets_from_cache {
            " (from index)"
        } else {
            ""
        },
        timings.scan_s,
        timings.assess_s,
        timings.table_s,
        t.elapsed().as_secs_f64()
    );
    let h = &table.header;
    eprintln!(
        "{}: {} records scanned, {} kept, {} targets in {:.1}s",
        h.sample.sample_id,
        h.records_scanned,
        h.records_kept,
        h.targets,
        started.elapsed().as_secs_f64()
    );
    for (state, n) in &h.states {
        eprintln!("  {state}: {n}");
    }
    Ok(())
}

fn is_cohort(path: &Path) -> Result<bool> {
    let mut magic = [0u8; 8];
    let mut f = std::fs::File::open(path).map_err(Error::io(path))?;
    let n = std::io::Read::read(&mut f, &mut magic).map_err(Error::io(path))?;
    Ok(n == 8 && &magic == pgsum::cohort::MAGIC)
}

fn extract_cohort(
    vcf: &Path,
    reference: &Path,
    packs: &[PathBuf],
    out: &Path,
    options: &pgsum::cohort::CohortOptions,
) -> Result<()> {
    let started = std::time::Instant::now();
    let reference = Reference::open(reference)?;
    let identity = reference_identity(&reference)?;
    let summary = pgsum::cohort::extract_cohort(vcf, &reference, &identity, packs, options, out)?;
    let h = &summary.header;
    eprintln!(
        "{} samples, {} records scanned, {} targets in {:.1}s",
        h.samples.len(),
        h.records_scanned,
        h.targets,
        started.elapsed().as_secs_f64()
    );
    for (state, n) in &h.states {
        eprintln!("  {state}: {n}");
    }
    Ok(())
}

/// Score every pack for every sample of a cohort file: `cohort-scores.tsv.zst`, one line per sample and score.
fn score_cohort(cohort: &Path, packs: &[PathBuf], out: &Path, options: &pgsum::score::Options) -> Result<()> {
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let started = std::time::Instant::now();
    let table = pgsum::cohort::CohortTable::open(cohort)?;
    table.preload()?;
    let scores: Vec<Result<pgsum::cohort::CohortScore>> = packs
        .par_iter()
        .map(|p| pgsum::cohort::score_cohort(&Pack::open(p)?, &table, options))
        .collect();
    let mut scores = scores.into_iter().collect::<Result<Vec<_>>>()?;
    scores.sort_by(|a, b| a.pgs_id.cmp(&b.pgs_id));
    let path = out.join("cohort-scores.tsv.zst");
    let tmp = path.with_extension("zst.tmp");
    let written = (|| -> std::io::Result<()> {
        let mut w = zstd::Encoder::new(BufWriter::new(std::fs::File::create(&tmp)?), 3)?;
        writeln!(
            w,
            "sample\tpgs_id\tpartial_raw_score\tscorable_terms\ttotal_terms\tterm_coverage\tweight_coverage\t\
             meets_coverage_guideline"
        )?;
        for (s, sample) in table.header.samples.iter().enumerate() {
            for c in &scores {
                let (tf, wf) = (
                    if c.total_terms == 0 {
                        0.0
                    } else {
                        c.scorable[s] as f64 / c.total_terms as f64
                    },
                    if c.effect_all == 0.0 {
                        0.0
                    } else {
                        c.effect_scorable[s] / c.effect_all
                    },
                );
                let meets = tf >= pgsum::score::COVERAGE_GUIDELINE && wf >= pgsum::score::COVERAGE_GUIDELINE;
                writeln!(
                    w,
                    "{sample}\t{}\t{}\t{}\t{}\t{tf:.6}\t{wf:.6}\t{}",
                    c.pgs_id,
                    c.text(s),
                    c.scorable[s],
                    c.total_terms,
                    meets as u8
                )?;
            }
        }
        w.finish()?.flush()?;
        std::fs::rename(&tmp, &path)
    })();
    written.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        Error::io(&path)(e)
    })?;
    eprintln!(
        "scored {} packs for {} samples in {:.1}s",
        scores.len(),
        table.header.samples.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Above this many packs, `score` decompresses the whole genotype table up front.
const PRELOAD_PACKS: usize = 50;

/// Score every pack against the table: packs in parallel, and each pack's terms in parallel chunks.
/// A reference panel, its groups, and the sample's nearest group.
struct Panel {
    name: String,
    cohort: pgsum::cohort::CohortTable,
    groups: pgsum::panel::Groups,
    ancestry: pgsum::panel::Ancestry,
}

fn open_panel(args: &ReferenceArgs, table: &GenotypeTable) -> Result<Option<Panel>> {
    let (Some(path), Some(groups)) = (&args.reference_panel, &args.reference_groups) else {
        return Ok(None);
    };
    let started = std::time::Instant::now();
    let cohort = pgsum::cohort::CohortTable::open(path)?;
    cohort.preload()?;
    let groups = pgsum::panel::read_groups(groups, &cohort, &args.reference_group_column)?;
    table.preload()?;
    let ancestry = pgsum::panel::nearest_group(table, &cohort, &groups)?;
    eprintln!(
        "reference panel {}: {} samples; nearest {} {} ({} sites, {:.1}s)",
        path.display(),
        cohort.header.samples.len(),
        groups.column,
        ancestry.nearest_group,
        ancestry.sites,
        started.elapsed().as_secs_f64()
    );
    Ok(Some(Panel {
        name: path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        cohort,
        groups,
        ancestry,
    }))
}

#[allow(clippy::too_many_arguments)]
fn score(
    table: &GenotypeTable,
    packs: &[PathBuf],
    out: &Path,
    terms: bool,
    bundle: bool,
    options: &pgsum::score::Options,
    panel: Option<&Panel>,
) -> Result<()> {
    std::fs::create_dir_all(out).map_err(Error::io(out))?;
    let started = std::time::Instant::now();
    // Blocks load as scoring reaches them; with many packs nearly all are needed, so load them in parallel.
    if packs.len() > PRELOAD_PACKS {
        table.preload()?;
    }
    let outcomes: Vec<(&PathBuf, Result<pgsum::score::ScoreResult>)> = packs
        .par_iter()
        .map(|path| {
            let result = (|| -> Result<pgsum::score::ScoreResult> {
                let pack = Pack::open(path)?;
                let id = &pack.header.pgs_id;
                let result = if terms {
                    let tsv = out.join(format!("{id}.terms.tsv"));
                    let mut w = BufWriter::new(std::fs::File::create(&tsv).map_err(Error::io(&tsv))?);
                    let r = pgsum::score::score(&pack, table, options, Some(&mut w))?;
                    w.flush().map_err(Error::io(&tsv))?;
                    r
                } else {
                    pgsum::score::score(&pack, table, options, None)?
                };
                let mut result = result;
                if let Some(p) = panel {
                    result.reference = Some(pgsum::panel::place(
                        &pack,
                        table,
                        &p.cohort,
                        &p.name,
                        &p.groups,
                        Some(&p.ancestry),
                        options,
                    )?);
                    result.ancestry = Some(p.ancestry.clone());
                }
                if !bundle {
                    let json_path = out.join(format!("{id}.score.json"));
                    let mut json = serde_json::to_vec_pretty(&result).map_err(|e| Error::Invalid(e.to_string()))?;
                    json.push(b'\n');
                    std::fs::write(&json_path, json).map_err(Error::io(&json_path))?;
                }
                Ok(result)
            })();
            (path, result)
        })
        .collect();
    let mut results = Vec::new();
    let mut failures = Vec::new();
    for (path, r) in outcomes {
        match r {
            Ok(r) => results.push(r),
            Err(e) => failures.push(format!("{}: {e}", path.display())),
        }
    }
    results.sort_by(|a, b| a.pgs_id.cmp(&b.pgs_id));
    let mut summary = String::from(
        "pgs_id\tstatus\traw_score\tpartial_raw_score\tscorable_terms\ttotal_terms\tterm_coverage\tweight_coverage\t\
         meets_coverage_guideline\tchrx_terms\treference_group\treference_percentile\tmatched_term_coverage\t\
         matched_weight_coverage\treference_meets_coverage_guideline\n",
    );
    for r in &results {
        let c = &r.partial.coverage;
        summary.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.pgs_id,
            r.status,
            r.raw_score.as_deref().unwrap_or(""),
            r.partial.raw_score,
            c.scorable_terms,
            c.total_terms,
            c.term_fraction,
            c.weight_fraction,
            r.partial.meets_coverage_guideline as u8,
            r.sex_chromosomes.chrx_terms,
            r.reference
                .as_ref()
                .and_then(|p| p.nearest_group.clone())
                .unwrap_or_default(),
            r.reference
                .as_ref()
                .and_then(|p| p.nearest_group_percentile)
                .map(|v| format!("{v:.2}"))
                .unwrap_or_default(),
            r.reference
                .as_ref()
                .map(|p| format!("{:.6}", p.matched_term_fraction))
                .unwrap_or_default(),
            r.reference
                .as_ref()
                .map(|p| format!("{:.6}", p.matched_weight_fraction))
                .unwrap_or_default(),
            r.reference
                .as_ref()
                .map(|p| (p.meets_coverage_guideline as u8).to_string())
                .unwrap_or_default()
        ));
        eprintln!(
            "{}: {} (partial {} over {:.2}% of terms, {:.2}% of weight)",
            r.pgs_id,
            r.status,
            r.partial.raw_score,
            100.0 * c.term_fraction,
            100.0 * c.weight_fraction
        );
    }
    if bundle {
        write_bundle(&out.join("results.jsonl.zst"), &results)?;
    }
    let summary_path = out.join("scores.tsv");
    std::fs::write(&summary_path, summary).map_err(Error::io(&summary_path))?;
    let terms_scored: u64 = results.iter().map(|r| r.required_terms).sum();
    eprintln!(
        "scored {} packs, {} terms, in {:.1}s",
        results.len(),
        terms_scored,
        started.elapsed().as_secs_f64()
    );
    if failures.is_empty() {
        Ok(())
    } else {
        Err(Error::Invalid(failures.join("\n")))
    }
}

/// Every result as one compact JSON line, zstd-compressed, written to a temporary file and renamed into place.
fn write_bundle(path: &Path, results: &[pgsum::score::ScoreResult]) -> Result<()> {
    let tmp = path.with_extension("zst.tmp");
    let written = (|| -> std::io::Result<()> {
        let mut w = zstd::Encoder::new(BufWriter::new(std::fs::File::create(&tmp)?), 3)?;
        for r in results {
            serde_json::to_writer(&mut w, r)?;
            w.write_all(b"\n")?;
        }
        w.finish()?.flush()?;
        std::fs::rename(&tmp, path)
    })();
    written.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        Error::io(path)(e)
    })
}

/// One line per target of a genotype table: position, alleles (a sequence target's normalized alleles; `*` for
/// a target counting any non-reference allele), call state and ALT dosage.
fn write_targets_tsv(table: &GenotypeTable, out: &mut impl Write) -> Result<()> {
    let io = |e| Error::io("<stdout>")(e);
    writeln!(out, "contig\tpos\tref\talt\tstate\talt_dosage\trefcall_adapted\tphased").map_err(io)?;
    for e in table.entries()? {
        let (code, pos) = pgsum::genotypes::key_position(e.key);
        let contig = pgsum::term::CONTIGS[code as usize - 1];
        let (pos, r, a) = match table.sequence(e.key) {
            Some(v) => (v.pos, v.ref_allele.clone(), v.alt.clone()),
            None => {
                let (_, _, r, a) = pgsum::genotypes::unpack_key(e.key);
                (pos as u64, (r as char).to_string(), (a as char).to_string())
            }
        };
        writeln!(
            out,
            "{contig}\t{pos}\t{r}\t{a}\t{}\t{}\t{}\t{}",
            e.state.as_str(),
            e.alt_dosage.map_or(String::new(), |d| d.to_string()),
            e.refcall_adapted as u8,
            e.phased as u8
        )
        .map_err(io)?;
    }
    Ok(())
}

fn evidence(pack: &Path, genotypes: &Path, reference: &Path, policy: Option<&str>, after: usize) -> Result<()> {
    let (pack, table, reference) = (
        Pack::open(pack)?,
        GenotypeTable::open(genotypes)?,
        Reference::open(reference)?,
    );
    let policy = policy.map_or_else(|| table.header.policy.id.clone(), str::to_owned);
    let mut out = BufWriter::new(std::io::stdout().lock());
    pgsum::evidence::write_rows(&pack, &table, &reference, &policy, after, &mut out)?;
    out.flush().map_err(Error::io("<stdout>"))?;
    Ok(())
}

fn inspect(path: &Path, header: bool, genotypes: Option<&Path>, options: &pgsum::score::Options) -> Result<()> {
    let mut out = BufWriter::new(std::io::stdout().lock());
    let json = |out: &mut BufWriter<_>, value: &dyn erased::Json| value.write(out);
    let is_table = std::fs::File::open(path)
        .and_then(|mut f| {
            let mut magic = [0u8; 8];
            std::io::Read::read_exact(&mut f, &mut magic).map(|_| &magic == pgsum::genotypes::MAGIC)
        })
        .map_err(Error::io(path))?;
    if is_table {
        let table = GenotypeTable::open(path)?;
        if header {
            json(&mut out, &table.header)?;
        } else {
            write_targets_tsv(&table, &mut out)?;
        }
    } else {
        let pack = Pack::open(path)?;
        if header {
            json(&mut out, &pack.header)?;
        } else {
            match genotypes {
                Some(g) => {
                    pgsum::score::score(&pack, &GenotypeTable::open(g)?, options, Some(&mut out))?;
                }
                None => pack.write_terms_tsv(&mut out)?,
            }
        }
    }
    out.flush().map_err(Error::io("<stdout>"))
}

mod erased {
    use pgsum::{Error, Result};

    pub trait Json {
        fn write(&self, out: &mut dyn std::io::Write) -> Result<()>;
    }

    impl<T: serde::Serialize> Json for T {
        fn write(&self, out: &mut dyn std::io::Write) -> Result<()> {
            serde_json::to_writer_pretty(&mut *out, self).map_err(|e| Error::Invalid(e.to_string()))?;
            writeln!(out).map_err(Error::io("<stdout>"))
        }
    }
}
