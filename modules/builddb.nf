// Build the shared search databases once; every mining channel reuses them from S3.

process FOLDSEEK_DB {
    label 'xlarge'
    conda 'bioconda::foldseek=9.427df8a'
    publishDir "${params.dbdir}/foldseek", mode: params.publish_mode

    input:
    path pdbs

    output:
    path 'humanAF*', emit: db

    script:
    """
    mkdir -p indir
    for f in ${pdbs}; do cp \$f indir/; done
    foldseek createdb indir humanAF --threads ${task.cpus}
    foldseek createindex humanAF tmp --threads ${task.cpus} || echo "createindex skipped"
    rm -rf tmp indir
    ls -la humanAF* | head
    """
}

process MMSEQS_DB {
    label 'large'
    conda 'bioconda::mmseqs2=15.6f452'
    publishDir "${params.dbdir}/mmseqs", mode: params.publish_mode

    input:
    path fasta

    output:
    path 'humanSeq*', emit: db

    script:
    """
    mmseqs createdb ${fasta} humanSeq
    ls -la humanSeq*
    """
}

process PLDDT_SUMMARY {
    label 'large'
    conda 'conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy'
    publishDir "${params.dbdir}", mode: params.publish_mode

    input:
    path pdbs

    output:
    path 'plddt_summary.csv', emit: csv

    script:
    """
    #!/usr/bin/env python
    import gzip, glob, csv, statistics, re, os

    rows = []
    for fp in sorted(glob.glob('*.pdb.gz')):
        m = re.search(r'AF-([A-Z0-9]+)-F(\\d+)-model', os.path.basename(fp))
        if not m:
            continue
        acc, frag = m.group(1), int(m.group(2))
        vals = []
        op = gzip.open(fp, 'rt')
        for line in op:
            # one B-factor (=pLDDT) per residue, read at CA only
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
        op.close()
        if not vals:
            continue
        rows.append({
            'accession': acc, 'fragment': frag, 'n_residues': len(vals),
            'plddt_mean': round(statistics.fmean(vals), 2),
            'plddt_median': round(statistics.median(vals), 2),
            'frac_plddt_gt70': round(sum(v > 70 for v in vals) / len(vals), 4),
            'frac_plddt_gt90': round(sum(v > 90 for v in vals) / len(vals), 4),
        })

    with open('plddt_summary.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('structures summarised:', len(rows))
    """
}
