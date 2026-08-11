#!/usr/bin/env nextflow
// phemine-nf — entry points for each stage of the Phe-binder mining campaign.
// Run one stage at a time with -entry, so a failure never forces a full restart.

include { FETCH_PROTEOME_FASTA; FETCH_PROTEOME_ANNOTATION;
          FETCH_AFDB; FETCH_PFAM }              from './modules/stage.nf'
include { FOLDSEEK_DB; MMSEQS_DB; PLDDT_SUMMARY } from './modules/builddb.nf'

// Stage 1 — pull reference data into S3 from inside AWS.
workflow stage_data {
    FETCH_PROTEOME_FASTA()
    FETCH_PROTEOME_ANNOTATION()
    FETCH_AFDB()
    FETCH_PFAM()
}

// Stage 2 — build the shared search databases from the staged data.
workflow build_db {
    fasta = Channel.fromPath("${params.outdir}/human_proteome.fasta", checkIfExists: true)
    pdbs  = Channel.fromPath("${params.outdir}/afdb/afdb_pdb/*.pdb.gz", checkIfExists: true)

    MMSEQS_DB(fasta)
    FOLDSEEK_DB(pdbs.collect())
    PLDDT_SUMMARY(pdbs.collect())
}

workflow {
    error "Specify a stage with -entry, e.g. -entry stage_data or -entry build_db"
}
