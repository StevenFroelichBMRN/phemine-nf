// Data staging: fetch reference data directly into S3 from inside AWS.
// Nothing large transits the orchestrating sandbox.

process FETCH_PROTEOME_FASTA {
    label 'network'
    conda 'conda-forge::curl conda-forge::gzip'
    publishDir "${params.outdir}", mode: params.publish_mode

    output:
    path 'human_proteome.fasta', emit: fasta

    script:
    def q = params.uniprot_query.replaceAll(' ', '+')
    """
    curl -sSfL --retry 5 --retry-delay 10 \\
      "https://rest.uniprot.org/uniprotkb/stream?query=${q}&format=fasta&compressed=true" \\
      -o proteome.fasta.gz
    gunzip -c proteome.fasta.gz > human_proteome.fasta
    n=\$(grep -c '^>' human_proteome.fasta)
    echo "sequences: \$n"
    # A truncated stream is the main failure mode here — fail loudly rather than
    # silently mining a partial proteome.
    if [ "\$n" -lt 15000 ]; then
        echo "ERROR: only \$n sequences; expected >15000 for reviewed human" >&2
        exit 1
    fi
    """
}

process FETCH_PROTEOME_ANNOTATION {
    label 'network'
    conda 'conda-forge::curl conda-forge::gzip'
    publishDir "${params.outdir}", mode: params.publish_mode

    output:
    path 'human_proteome_annotation.tsv', emit: tsv

    script:
    def q = params.uniprot_query.replaceAll(' ', '+')
    def fields = [
        'accession','id','protein_name','gene_primary','length','ec',
        // NB: 'ft_metal' was retired from the UniProt API (400 Invalid fields parameter);
        // metal-binding is reported under ft_binding now.
        'cc_catalytic_activity','cc_cofactor','ft_binding','ft_act_site','ft_site',
        'ft_transmem','cc_subcellular_location','go_f','go_p',
        'xref_pfam','xref_interpro','xref_supfam','xref_gene3d',
        'xref_pdb','xref_alphafolddb','xref_chembl','xref_drugbank',
        'cc_similarity','protein_families','rhea','cc_function'
    ].join(',')
    """
    curl -sSfL --retry 5 --retry-delay 10 \\
      "https://rest.uniprot.org/uniprotkb/stream?query=${q}&format=tsv&compressed=true&fields=${fields}" \\
      -o ann.tsv.gz
    gunzip -c ann.tsv.gz > human_proteome_annotation.tsv
    n=\$(( \$(wc -l < human_proteome_annotation.tsv) - 1 ))
    echo "annotation rows: \$n"
    if [ "\$n" -lt 15000 ]; then
        echo "ERROR: only \$n annotation rows" >&2
        exit 1
    fi
    """
}

process FETCH_AFDB {
    label 'medium'
    // 5.1 GB tar -> ~23k gzipped PDBs; the default Batch scratch is not enough.
    disk 120.GB
    conda 'conda-forge::curl conda-forge::tar conda-forge::gzip'
    publishDir "${params.outdir}/afdb", mode: params.publish_mode

    output:
    path 'afdb_pdb/*.pdb.gz', emit: structures
    path 'afdb_manifest.txt',  emit: manifest

    script:
    """
    curl -sSfL --retry 5 --retry-delay 15 "${params.afdb_url}" -o afdb.tar
    mkdir -p afdb_raw afdb_pdb
    tar -xf afdb.tar -C afdb_raw
    # Keep PDB models only; drop the PAE JSONs and CIFs to cut S3 footprint.
    find afdb_raw -name '*model_v4.pdb.gz' -exec mv {} afdb_pdb/ \\;
    ls afdb_pdb | wc -l > afdb_manifest.txt
    ls afdb_pdb >> afdb_manifest.txt
    n=\$(ls afdb_pdb | wc -l)
    echo "AFDB models: \$n"
    if [ "\$n" -lt 15000 ]; then
        echo "ERROR: only \$n AFDB models extracted" >&2
        exit 1
    fi
    rm -f afdb.tar
    rm -rf afdb_raw
    """
}

process FETCH_PFAM {
    label 'network'
    conda 'conda-forge::curl conda-forge::gzip'
    publishDir "${params.outdir}/pfam", mode: params.publish_mode

    output:
    path 'Pfam-A.hmm', emit: hmm

    script:
    """
    curl -sSfL --retry 5 --retry-delay 15 "${params.pfam_url}" -o Pfam-A.hmm.gz
    gunzip -c Pfam-A.hmm.gz > Pfam-A.hmm
    grep -c '^NAME' Pfam-A.hmm
    """
}
