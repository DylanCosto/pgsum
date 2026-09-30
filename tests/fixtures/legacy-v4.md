`legacy-v4.pgsg` was generated with the pre-expansion pgsum v0.5.1 executable
(baseline from commit 1f419871a70446b5fbeb68ae83e4a56d4d81ae95), using the existing
synthetic.fa, synthetic.g.vcf.gz, and PGS999998_hmPOS_GRCh38.txt.gz plus metadata.
It contains only synthetic data. Reproduction with that executable:

```
pgsum compile PGS999998_hmPOS_GRCh38.txt.gz --reference synthetic.fa --out packs/
pgsum extract --gvcf synthetic.g.vcf.gz --reference synthetic.fa --pack packs/ --out legacy-v4.pgsg
```

The compatibility test reads the checked-in v4 file, rewrites v5, and verifies every
call and source record. Do not regenerate this fixture with the current writer.
