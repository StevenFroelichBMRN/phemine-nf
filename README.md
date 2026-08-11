# phemine-nf

Nextflow pipeline for a proteome-wide search for human proteins that could be
engineered to bind and transform L-phenylalanine.

Heavy compute runs on AWS Batch via Seqera Platform; the orchestrating session
holds only compact result tables.

## Stages

Run one stage at a time with `-entry`:

| entry | what it does | scale |
|---|---|---|
| `stage_data` | fetch UniProt human proteome, AFDB structures, Pfam-A into S3 | network-bound |
| `build_db` | build Foldseek + MMseqs2 databases and a pLDDT summary | 32 cpu |

## Layout

All data lives under `s3://<bucket>/phemine/`:

```
raw/    proteome fasta, annotation tsv, afdb/, pfam/
db/     foldseek/, mmseqs/, plddt_summary.csv
channels/  per-channel hit tables and scores
pockets/   fpocket features
```

Each stage publishes to its own prefix, so a re-run skips work whose output
already exists.
