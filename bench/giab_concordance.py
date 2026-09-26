#!/usr/bin/env python3
"""Compare pgsum's calls at every target of a genotype table with a GIAB truth set.

Usage:
    pgsum inspect sample.pgsg | python3 giab_concordance.py TRUTH.vcf.gz TRUTH.bed OUT_DIR

Only targets on chr1-chr22 inside the truth set's confident regions (BED) are evaluated. There, a site with
no truth record is homozygous reference. The truth dosage of a target is the number of copies of its ALT
(for an any-allele target `*`, of any non-reference allele) in the truth genotype. Targets where a truth
record other than a plain SNV at the same position touches the target's span (an indel at or across it, or
an MNP) are counted as `complex` and not compared, since their dosage depends on representation.

Writes OUT_DIR/summary.json and OUT_DIR/discordant.tsv (every compared target where the call and the truth
disagree).
"""

import bisect
import collections
import gzip
import json
import os
import sys

CHROMS = {f"chr{i}" for i in range(1, 23)}


def read_bed(path):
    starts, ends = collections.defaultdict(list), collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            c, s, e = line.split("\t")[:3]
            if c in CHROMS:
                starts[c].append(int(s))
                ends[c].append(int(e))
    for c in starts:
        order = sorted(range(len(starts[c])), key=starts[c].__getitem__)
        starts[c] = [starts[c][i] for i in order]
        ends[c] = [ends[c][i] for i in order]
    return starts, ends


def in_bed(bed, c, pos):
    """1-based `pos` inside a 0-based half-open interval."""
    starts, ends = bed
    i = bisect.bisect_right(starts[c], pos - 1) - 1
    return i >= 0 and ends[c][i] >= pos


def read_truth(path):
    """Per chromosome: records by POS, and the positions covered by a record beyond its POS."""
    records = collections.defaultdict(dict)
    spanned = collections.defaultdict(set)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alt, _, _, _, fmt, sample = line.rstrip("\n").split("\t")[:10]
            if c not in CHROMS:
                continue
            pos = int(pos)
            gt = sample.split(":")[fmt.split(":").index("GT")]
            gt = tuple(int(a) for a in gt.replace("|", "/").split("/") if a != ".")
            alts = tuple(alt.split(","))
            records[c].setdefault(pos, []).append((ref, alts, gt))
            if len(ref) > 1:
                spanned[c].update(range(pos + 1, pos + len(ref)))
    return records, spanned


def truth_dosage(records, spanned, c, pos, ref, alt):
    """Truth copies of `alt` for the target, 'complex', or 'ref_mismatch'."""
    span = range(pos, pos + len(ref))
    if any(p in spanned[c] for p in span):
        return "complex"
    here = records[c].get(pos, [])
    if any(records[c].get(p) for p in span[1:]):
        return "complex"
    if not here:
        return 0
    snv = alt == "*" or (len(ref) == 1 and len(alt) == 1)
    dosage = 0
    for tref, talts, gt in here:
        alleles = (tref,) + talts
        plain_snv = len(tref) == 1 and all(len(a) == 1 for a in talts)
        if snv:
            if not plain_snv:
                return "complex"
            if tref != ref:
                return "ref_mismatch"
            if alt == "*":
                dosage += sum(1 for a in gt if a != 0)
            else:
                dosage += sum(1 for a in gt if alleles[a] == alt)
        else:
            if tref == ref and alt in talts:
                dosage += sum(1 for a in gt if alleles[a] == alt)
            else:
                return "complex"
    return min(dosage, 2)


def main():
    truth_vcf, truth_bed, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    bed = read_bed(truth_bed)
    records, spanned = read_truth(truth_vcf)
    counts = collections.Counter()
    matrix = {kind: collections.Counter() for kind in ("snv", "sequence", "any_allele")}
    uncalled = collections.Counter()
    other_allele = collections.Counter()
    discordant = open(os.path.join(out_dir, "discordant.tsv"), "w")
    discordant.write("contig\tpos\tref\talt\tstate\tpgsum_dosage\ttruth_dosage\n")
    header = sys.stdin.readline().rstrip("\n").split("\t")
    assert header[:6] == ["contig", "pos", "ref", "alt", "state", "alt_dosage"], header
    for line in sys.stdin:
        c, pos, ref, alt, state, dosage = line.rstrip("\n").split("\t")[:6]
        counts["targets"] += 1
        if c not in CHROMS:
            continue
        counts["targets_chr1_22"] += 1
        pos = int(pos)
        if not in_bed(bed, c, pos):
            continue
        counts["in_confident_regions"] += 1
        truth = truth_dosage(records, spanned, c, pos, ref, alt)
        if not isinstance(truth, int):
            counts[f"truth_{truth}"] += 1
            continue
        counts["comparable"] += 1
        kind = "any_allele" if alt == "*" else "snv" if len(ref) == 1 and len(alt) == 1 else "sequence"
        if state == "other_called_allele":
            other_allele["truth_has_target_alt" if truth else "truth_lacks_target_alt"] += 1
        if dosage == "":
            uncalled[state] += 1
            continue
        called = int(dosage)
        matrix[kind][f"truth{truth}_pgsum{called}"] += 1
        if called != truth:
            discordant.write(f"{c}\t{pos}\t{ref}\t{alt}\t{state}\t{called}\t{truth}\n")
    discordant.close()
    summary = {"counts": dict(counts), "uncalled_by_state": dict(uncalled), "other_called_allele": dict(other_allele)}
    for kind, m in matrix.items():
        total = sum(m.values())
        agree = sum(v for k, v in m.items() if k[5] == k[-1])
        nonref = sum(v for k, v in m.items() if not (k[5] == "0" and k[-1] == "0"))
        nonref_agree = sum(v for k, v in m.items() if k[5] == k[-1] != "0")
        summary[kind] = {
            "called": total,
            "concordant": agree,
            "concordance": agree / total if total else None,
            "non_reference": nonref,
            "non_reference_concordance": nonref_agree / nonref if nonref else None,
            "matrix": dict(sorted(m.items())),
        }
    with open(os.path.join(out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
