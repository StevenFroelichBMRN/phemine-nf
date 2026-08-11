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

include { PHMMER_SEARCH; JACKHMMER_SEARCH;
          HMMSEARCH_PFAM; MMSEQS_PROFILE_SEARCH }  from './modules/homology.nf'
include { FOLDSEEK_SEARCH; FPOCKET_SCAN; SMINA_DOCK_PHE } from './modules/structure.nf'

// Stage 3 — sequence remote-homology channel.
workflow homology {
    proteome = Channel.fromPath("${params.outdir}/human_proteome.fasta", checkIfExists: true)
    queries  = Channel.fromPath("${params.ref_fasta}/*.fasta", checkIfExists: true)

    PHMMER_SEARCH(queries.combine(proteome))
    JACKHMMER_SEARCH(queries.combine(proteome))
}

workflow pfam_scan {
    proteome = Channel.fromPath("${params.outdir}/human_proteome.fasta", checkIfExists: true)
    hmm      = Channel.fromPath("${params.outdir}/pfam/Pfam-A.hmm", checkIfExists: true)
    HMMSEARCH_PFAM(hmm.combine(proteome))
}

// Stage 4 — structural remote-homology channel.
workflow structure_search {
    queries = Channel.fromPath("${params.ref_structures}/*.pdb", checkIfExists: true).collect()
    db      = Channel.fromPath("${params.dbdir}/foldseek/humanAF*", checkIfExists: true).collect()
    FOLDSEEK_SEARCH(queries.combine(db))
}

// Stage 5 — pocket detection over the shortlist (one task per structure).
workflow pockets {
    structures = Channel.fromPath("${params.shortlist}", checkIfExists: true)
                        .splitCsv(header: true)
                        .map { row -> file("${params.outdir}/afdb/afdb_pdb/AF-${row.accession}-F1-model_v4.pdb.gz") }
                        .filter { it.exists() }
    FPOCKET_SCAN(structures)
}
