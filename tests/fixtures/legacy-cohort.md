# Legacy cohort reader fixtures

These synthetic files were written by pgsum commit `1521f9f`, before cohort diagnostics were added.
The writer binary SHA-256 was `89f6d3367a1629f068aaf6cd91af107cde93c1d5ae575e930c5847db0cd82528`.
Do not regenerate them with the current writer.

- `legacy-cohort-v1.pgsc`: SHA-256 `5914a1a6caf8e3d49b0a393e558536ab7dfe0857611565bb65869b2e79bef48d`.
- `legacy-cohort-v2.pgsc`: SHA-256 `b1e3e623332897a32b2b736c52930287d75c1b97c61a7135e318eb9eeb44216d`.

Compile the existing `PGS999997_hmPOS_GRCh38.txt.gz` and adjacent metadata against `synthetic.fa`.
The input `legacy.vcf` is the following tab-separated text:

```text
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
1	2	.	C	T	.	PASS	.	GT:DS:GP	0/1:0.4:0.25,0.5,0.25	./.:.:.	1/1:2:0,0,1
```

```sh
pgsum compile PGS999997_hmPOS_GRCh38.txt.gz --reference synthetic.fa --out packs/
pgsum extract --gvcf legacy.vcf --reference synthetic.fa --pack packs/ --all-samples --accept-missing-quality --threads 2 --out legacy-cohort-v1.pgsc
pgsum extract --gvcf legacy.vcf --reference synthetic.fa --pack packs/ --all-samples --accept-missing-quality --threads 2 --dosage-field gp --out legacy-cohort-v2.pgsc
```

Headers contain original paths, so reproduction in another directory need not have the same file hash.
The compatibility test checks preserved scores and explicitly unavailable legacy exclusion reasons.
