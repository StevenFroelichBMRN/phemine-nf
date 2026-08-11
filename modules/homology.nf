// Sequence-based remote homology: one task per reference query, fanned out on Batch.

process PHMMER_SEARCH {
    tag "${query.baseName}"
    label 'small'
    conda 'bioconda::hmmer=3.4'
    publishDir "${params.s3_base}/channels/homology/phmmer", mode: params.publish_mode

    input:
    tuple path(query), path(proteome)

    output:
    path "${query.baseName}.phmmer.tbl", emit: tbl

    script:
    """
    phmmer --cpu ${task.cpus} -E 10 --domE 10 \\
           --tblout ${query.baseName}.phmmer.tbl \\
           --domtblout ${query.baseName}.phmmer.dom \\
           ${query} ${proteome} > /dev/null
    """
}

process JACKHMMER_SEARCH {
    tag "${query.baseName}"
    label 'medium'
    conda 'bioconda::hmmer=3.4'
    publishDir "${params.s3_base}/channels/homology/jackhmmer", mode: params.publish_mode

    input:
    tuple path(query), path(proteome)

    output:
    path "${query.baseName}.jack.tbl", emit: tbl

    script:
    """
    # Iterated search recovers homologs too diverged for a single-pass phmmer.
    jackhmmer --cpu ${task.cpus} -N 3 -E 10 --domE 10 \\
              --tblout ${query.baseName}.jack.tbl \\
              ${query} ${proteome} > /dev/null
    """
}

process HMMSEARCH_PFAM {
    label 'xlarge'
    conda 'bioconda::hmmer=3.4'
    publishDir "${params.s3_base}/channels/homology/pfam", mode: params.publish_mode

    input:
    tuple path(hmm), path(proteome)

    output:
    path 'pfam_domains.tbl', emit: tbl

    script:
    """
    hmmsearch --cpu ${task.cpus} --cut_ga \\
              --domtblout pfam_domains.tbl \\
              ${hmm} ${proteome} > /dev/null
    """
}

process MMSEQS_PROFILE_SEARCH {
    label 'xlarge'
    conda 'bioconda::mmseqs2=15.6f452'
    publishDir "${params.s3_base}/channels/homology/mmseqs", mode: params.publish_mode

    input:
    tuple path(ref_fasta), path(dbfiles)

    output:
    path 'mmseqs_hits.m8', emit: hits

    script:
    """
    mmseqs createdb ${ref_fasta} refDB
    # Profile search: build profiles from the reference queries, then search the
    # human DB with them — more sensitive than sequence-sequence at low identity.
    mmseqs search refDB humanSeq resDB tmp \\
        --threads ${task.cpus} -s 7.5 --max-seqs 4000 -e 10 --num-iterations 3
    mmseqs convertalis refDB humanSeq resDB mmseqs_hits.m8 --threads ${task.cpus} \\
        --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'
    rm -rf tmp
    """
}
