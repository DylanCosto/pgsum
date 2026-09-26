"""Write the synthetic genotype fixture: PGS999998 (one SNV term per gVCF scenario) and synthetic.g.vcf.gz.

Run from this directory; needs bgzip and tabix. The reference is synthetic.fa (60 bp chr1, 30 bp chrX).
"""
import gzip
import json
import subprocess

chr1 = "ACGTTGCAAC" + "NNacgttgca" + "ACGTTGCAAC" * 4


def seq(p):
    return chr1[p - 1].upper()


# (pos, effect, other, weight); every term resolves on the forward strand.
terms = [(1, "G", "A", "0.1"), (2, "C", "T", "0.2"), (3, "G", "A", "-0.3"), (4, "C", "T", "0.4"), (5, "G", "T", "0.5"),
         (6, "A", "G", "0.6"), (7, "T", "C", "0.7"), (8, "G", "A", "0.8"), (9, "C", "A", "0.9"), (10, "T", "C", "1.0"),
         (13, "G", "A", "1.1"), (14, "A", "C", "1.2"), (15, "A", "G", "1.3"), (16, "G", "T", "1.4"), (17, "C", "T", "1.5"),
         (18, "A", "G", "1.6"), (19, "T", "C", "1.7"), (21, "G", "A", "1.8"), (22, "T", "C", "1.9"), (23, "A", "G", "2.0"),
         (24, "G", "T", "2.1"), (25, "G", "T", "2.2"), (26, "A", "G", "2.3"), (27, "T", "C", "2.4"), (28, "G", "A", "2.5")]
for p, e, o, _ in terms:
    assert seq(p) in (e, o) and {e, o} not in ({"A", "T"}, {"C", "G"}), (p, e, o, seq(p))
cols = ["rsID", "chr_name", "effect_allele", "other_allele", "effect_weight", "hm_source", "hm_chr", "hm_pos"]
rows = ["\t".join([f"rs{p}", "1", e, o, w, "ENSEMBL", "1", str(p)]) for p, e, o, w in terms]
rows.append("\t".join(["rsX2", "X", "C", "A", "3.0", "ENSEMBL", "X", "2"]))
header = ["###PGS CATALOG SCORING FILE - synthetic test fixture", "#format_version=2.0", "#pgs_id=PGS999998",
          f"#variants_number={len(rows)}", "#weight_type=beta", "#HmPOS_build=GRCh38"]
with gzip.GzipFile("PGS999998_hmPOS_GRCh38.txt.gz", "wb", compresslevel=9, mtime=0) as f:
    f.write(("\n".join(header + ["\t".join(cols)] + rows) + "\n").encode())
with open("PGS999998.metadata.json", "w") as f:
    json.dump(dict(id="PGS999998", name="synthetic genotypes", variants_number=len(rows), weight_type="beta",
                   matches_publication=True, license="CC0 1.0 (synthetic test fixture)"), f, indent=1)


def block(pos, end, gt="0/0", gq="50", min_dp="30", ref=None, fmt="GT:GQ:MIN_DP:PL", sample=None):
    return ["chr1", pos, ref or seq(pos), "<*>", "0", ".", f"END={end}", fmt, sample or f"{gt}:{gq}:{min_dp}:0,60,900"]


def variant(pos, ref, alt, gt, filt="PASS", gq="40", dp="30", fmt="GT:GQ:DP:AD:VAF:PL", sample=None, chrom="chr1"):
    return [chrom, pos, ref, alt, "30.5", filt, ".", fmt, sample or f"{gt}:{gq}:{dp}:15,15:0.5:30,0,40"]


haploid = "1:40:30:0,30:1:40,0"
records = [
    block(1, 3),                                                   # 1-3: confident reference block
    variant(4, "T", "C,<*>", "0/1"),                               # 4: heterozygous
    variant(5, "T", "G,<*>", "1/1"),                               # 5: homozygous ALT
    variant(6, "G", "A,<*>", "0/0", filt="RefCall"),               # 6: DeepVariant RefCall
    variant(7, "C", "T,<*>", "0/1", fmt="GT:GQ:DP:FT", sample="0/1:40:30:LowGQ"),  # 7: per-sample FT
    variant(8, "A", "G,<*>", "0/1", gq="15"),                      # 8: low GQ
    variant(9, "A", "C,<*>", "./."),                               # 9: no call
    variant(10, "C", "T,<*>", "./1"),                              # 10: partial no call
    variant(13, "A", "G,<*>", "1", sample=haploid),                # 13: haploid
    variant(14, "C", "A,<*>", "0/1", filt="LowQual"),              # 14: site filter
    block(15, 15, fmt="GT:GQ", sample="0/0:50"),                   # 15: block without MIN_DP
    variant(16, "T", "C,G,<*>", "0/1"),                            # 16: called ALT is not the term's ALT
    variant(17, "T", "G,C,<*>", "0/2"),                            # 17: term ALT is the second ALT
    variant(18, "GC", "G,<*>", "0/1"),                             # 18-19: deletion
    variant(20, "AA", "A,<*>", "0/1"),                             # 20-21: deletion ...
    block(21, 22),                                                 # ... and a block overlapping it at 21
    block(23, 23, ref="T"),                                        # 23: block REF differs from the reference
    variant(24, "T", "G,<*>", "0|1", fmt="GT:GQ:DP:PS", sample="0|1:40:30:1234"),  # 24: phased
    # 25: no record
    block(26, 26, gt="0/1"),                                       # 26: block with a non-reference genotype
    block(27, 27, gq="."),                                         # 27: block with GQ missing
    variant(28, "A", "G,<*>", "0/1", filt="."),                    # 28: FILTER "."
    variant(5, "A", "T,<*>", "0/1", chrom="chrUn_test"),           # other contigs are skipped
    variant(2, "A", "C,<*>", "1", chrom="chrX", sample=haploid),   # chrX haploid call
]
vcf_header = [
    "##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">',
    '##FILTER=<ID=LowQual,Description="Low quality">', '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in block">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">', "##DeepVariant_version=1.10.0",
    "##contig=<ID=chr1,length=60>", "##contig=<ID=chrUn_test,length=100>", "##contig=<ID=chrX,length=30>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH",
]
body = ["\t".join([r[0], str(r[1]), "."] + [str(x) for x in r[2:]]) for r in records]
with open("synthetic.g.vcf", "w") as f:
    f.write("\n".join(vcf_header + body) + "\n")
subprocess.run(["bgzip", "-f", "synthetic.g.vcf"], check=True)
subprocess.run(["tabix", "-f", "-p", "vcf", "synthetic.g.vcf.gz"], check=True)

# PGS999997: every term sits on a passing call above, with each weight model and a spread of exponents, so the
# complete-score path and exact sums are exercised.
complete = [
    # pos, effect, other, model, weights
    (1, "G", "A", "additive", ["0.1"]),
    (2, "C", "T", "dominant", ["-2.5E-3"]),
    (3, "G", "A", "recessive", ["1234"]),
    (4, "C", "T", "dosage", ["0", "0.15", "0.3"]),
    (5, "G", "T", "additive", ["0.000001"]),
    (6, "A", "G", "additive", ["-0.75"]),
    (17, "C", "T", "additive", ["1E+2"]),
    (22, "T", "C", "additive", ["0.12345678901234567890123"]),
    (24, "G", "T", "recessive", ["5"]),
    (28, "A", "G", "dominant", ["-1.5"]),
]
cols = ["rsID", "chr_name", "effect_allele", "other_allele", "effect_weight", "dosage_0_weight", "dosage_1_weight",
        "dosage_2_weight", "is_dominant", "is_recessive", "hm_source", "hm_chr", "hm_pos"]
rows = []
for p, e, o, model, w in complete:
    assert seq(p) in (e, o) and {e, o} not in ({"A", "T"}, {"C", "G"}), (p, e, o, seq(p))
    dosage = w if model == "dosage" else ["", "", ""]
    rows.append("\t".join([f"rs{p}", "1", e, o, "" if model == "dosage" else w[0], *dosage,
                           "TRUE" if model == "dominant" else "", "TRUE" if model == "recessive" else "",
                           "ENSEMBL", "1", str(p)]))
header = ["###PGS CATALOG SCORING FILE - synthetic test fixture", "#format_version=2.0", "#pgs_id=PGS999997",
          f"#variants_number={len(rows)}", "#weight_type=beta", "#HmPOS_build=GRCh38"]
with gzip.GzipFile("PGS999997_hmPOS_GRCh38.txt.gz", "wb", compresslevel=9, mtime=0) as f:
    f.write(("\n".join(header + ["\t".join(cols)] + rows) + "\n").encode())
with open("PGS999997.metadata.json", "w") as f:
    json.dump(dict(id="PGS999997", name="synthetic complete score", variants_number=len(rows), weight_type="beta",
                   matches_publication=True, license="CC0 1.0 (synthetic test fixture)"), f, indent=1)


# Scores without an author other allele, for the opt-in inference rules. Each row sits on a scenario above.
def effect_only(pgs_id, name, rows, infer_column=False):
    """rows: (pos, effect, weight[, hm_inferOtherAllele])."""
    cols = ["rsID", "chr_name", "effect_allele", "effect_weight", "hm_source", "hm_chr", "hm_pos"]
    if infer_column:
        cols.append("hm_inferOtherAllele")
    lines = []
    for r in rows:
        p, e, w = r[:3]
        line = [f"rs{p}", "1", e, w, "ENSEMBL", "1", str(p)]
        if infer_column:
            line.append(r[3])
        lines.append("\t".join(line))
    header = ["###PGS CATALOG SCORING FILE - synthetic test fixture", "#format_version=2.0", f"#pgs_id={pgs_id}",
              f"#variants_number={len(lines)}", "#weight_type=beta", "#HmPOS_build=GRCh38"]
    with gzip.GzipFile(f"{pgs_id}_hmPOS_GRCh38.txt.gz", "wb", compresslevel=9, mtime=0) as f:
        f.write(("\n".join(header + ["\t".join(cols)] + lines) + "\n").encode())
    with open(f"{pgs_id}.metadata.json", "w") as f:
        json.dump(dict(id=pgs_id, name=name, variants_number=len(lines), weight_type="beta", matches_publication=True,
                       license="CC0 1.0 (synthetic test fixture)"), f, indent=1)


# Every effect allele differs from the reference base: the other allele is the reference.
effect_only("PGS999996", "synthetic effect is ALT", [
    (1, "G", "0.1"), (4, "C", "0.2"), (5, "G", "0.3"), (16, "G", "0.4"), (17, "C", "0.5"), (22, "T", "0.6")])
# Every effect allele is the reference base: count reference copies (any non-reference allele counts against).
effect_only("PGS999995", "synthetic effect is REF", [
    (1, "A", "0.1"), (4, "T", "0.2"), (5, "T", "0.3"), (6, "G", "0.4"), (16, "T", "0.5"), (17, "T", "0.6")])
# No convention (half and half): only the Catalog's inferred allele can orient a term.
effect_only("PGS999994", "synthetic Catalog-inferred", [
    (4, "C", "0.2", "T"), (5, "T", "0.3", "G"), (2, "C", "0.4", "A/T"), (1, "A", "0.5", "T"), (24, "G", "0.6", "")],
    infer_column=True)


# Indels: a second gVCF over chr1:31-60 (ACGTTGCAAC ×3), a public variant set and a custom score.
def indel_block(pos, end):
    return ["chr1", pos, seq(pos), "<*>", "0", ".", f"END={end}", "GT:GQ:MIN_DP:PL", "0/0:50:30:0,60,900"]


indel_records = [
    indel_block(31, 32),
    variant(33, "GT", "G,<*>", "0/1"),          # 33-34: deletion of T, left-aligned
    indel_block(35, 42),
    variant(43, "GTT", "GT,<*>", "0/1"),        # 43-45: the same kind of deletion, written with extra context
    indel_block(46, 50),
    indel_block(51, 52),                        # 51-55: two adjacent blocks cover the 3-base deletion at 51
    indel_block(53, 55),
    indel_block(56, 60),
]
body = ["\t".join([r[0], str(r[1]), "."] + [str(x) for x in r[2:]]) for r in indel_records]
with open("synthetic_indel.g.vcf", "w") as f:
    f.write("\n".join(vcf_header + body) + "\n")
subprocess.run(["bgzip", "-f", "synthetic_indel.g.vcf"], check=True)
with open("public_variants.tsv", "w") as f:
    f.write("#synthetic public variant set\nchr1\t33\tGT\tG\nchr1\t43\tGT\tG\nchr1\t51\tACG\tA\n")
with open("INDELS.tsv", "w") as f:
    f.write("#pgs_id=INDELS\n#genome_build=GRCh38\n"
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
            "1\t33\tG\tGT\t0.5\n"      # public pair GT>G; effect is ALT; heterozygous
            "1\t43\tGT\tG\t0.25\n"     # public pair GT>G; effect is REF; heterozygous via GTT>GT
            "1\t51\tA\tACG\t2\n"       # public pair ACG>A; homozygous reference over two blocks
            "1\t56\tGA\tG\t1.5\n"      # only G fits the reference (GC): insertion G>GA; homozygous reference
            "1\t36\tG\tGC\t3\n")       # both fit (GC), no public record: stays unresolved
