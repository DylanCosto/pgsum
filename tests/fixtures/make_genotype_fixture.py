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
