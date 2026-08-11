// Structure-based channels: Foldseek remote-homology search and fpocket geometry.

process FOLDSEEK_SEARCH {
    label 'xlarge'
    conda 'bioconda::foldseek=9.427df8a'
    publishDir "${params.s3_base}/channels/structure", mode: params.publish_mode

    input:
    path queries, stageAs: 'qdir/*'
    path dbfiles

    output:
    path 'foldseek_hits.tsv', emit: hits

    script:
    """
    # Both inputs are collected lists; staging queries into qdir/ directly avoids
    # a combine() that would pair each query with the whole DB list.
    foldseek createdb qdir queryDB

    # Exhaustive-ish sensitivity: structural homologs of Phe enzymes are the point,
    # and they are frequently invisible to sequence search.
    foldseek search queryDB humanAF resDB tmp \\
        --threads ${task.cpus} -s 9.5 -e 10 --max-seqs 5000 -a

    foldseek convertalis queryDB humanAF resDB foldseek_hits.tsv \\
        --threads ${task.cpus} \\
        --format-output 'query,target,fident,alnlen,evalue,bits,qstart,qend,tstart,tend,qcov,tcov,alntmscore,qtmscore,ttmscore,lddt,prob,qaln,taln'
    rm -rf tmp
    wc -l foldseek_hits.tsv
    """
}

process FPOCKET_SCAN {
    tag "${pdb.simpleName}"
    label 'tiny'
    conda 'bioconda::fpocket=4.1'
    publishDir "${params.s3_base}/pockets/raw", mode: params.publish_mode

    input:
    path pdb

    output:
    path "${pdb.simpleName}_pockets.csv", optional: true, emit: csv

    script:
    """
    set +e
    if [[ "${pdb}" == *.gz ]]; then
        gunzip -c ${pdb} > model.pdb
    else
        cp ${pdb} model.pdb
    fi

    fpocket -f model.pdb > /dev/null 2>&1

    INFO=model_out/model_info.txt
    if [ ! -f "\$INFO" ]; then
        echo "no pockets for ${pdb.simpleName}" >&2
        exit 0
    fi

    python3 - <<'PY'
import re, csv, os
acc = "${pdb.simpleName}"
txt = open("model_out/model_info.txt").read()
blocks = re.split(r'Pocket\\s+(\\d+)\\s*:', txt)
rows = []
for i in range(1, len(blocks), 2):
    pid = int(blocks[i]); body = blocks[i+1]
    d = {}
    for line in body.strip().splitlines():
        if ':' in line:
            k, v = line.split(':', 1)
            try: d[k.strip()] = float(v.strip())
            except ValueError: pass
    rows.append({
        'structure': acc, 'pocket': pid,
        'score': d.get('Score'),
        'druggability': d.get('Druggability Score'),
        'n_alpha_spheres': d.get('Number of Alpha Spheres'),
        'volume': d.get('Volume'),
        'hydrophobicity': d.get('Hydrophobicity score'),
        'polarity': d.get('Polarity score'),
        'charge': d.get('Charge score'),
        'apolar_sasa_prop': d.get('Apolar alpha sphere proportion'),
        'mean_loc_hyd_dens': d.get('Mean local hydrophobic density'),
    })
if rows:
    with open(f"{acc}_pockets.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
PY
    """
}

process SMINA_DOCK_PHE {
    tag "${receptor.simpleName}"
    label 'small'
    conda 'bioconda::smina=2020.12.10 bioconda::openbabel=3.1.1'
    publishDir "${params.s3_base}/docking", mode: params.publish_mode

    input:
    tuple path(receptor), val(cx), val(cy), val(cz), path(ligand)

    output:
    path "${receptor.simpleName}_dock.csv", optional: true, emit: csv

    script:
    """
    set +e
    if [[ "${receptor}" == *.gz ]]; then gunzip -c ${receptor} > rec.pdb; else cp ${receptor} rec.pdb; fi
    obabel rec.pdb -O rec.pdbqt -xr 2>/dev/null

    # exhaustiveness 16: rank correlation with 64 measured at 1.000 on a prior probe,
    # so the extra cost buys nothing here.
    smina --receptor rec.pdbqt --ligand ${ligand} \\
          --center_x ${cx} --center_y ${cy} --center_z ${cz} \\
          --size_x 20 --size_y 20 --size_z 20 \\
          --exhaustiveness 16 --num_modes 5 --seed 42 --cpu ${task.cpus} \\
          --out poses.sdf --log dock.log 2>/dev/null

    python3 - <<'PY'
import re, csv
acc = "${receptor.simpleName}"
best = None
try:
    for line in open("dock.log"):
        m = re.match(r'\\s*(\\d+)\\s+(-?\\d+\\.\\d+)', line)
        if m and int(m.group(1)) == 1:
            best = float(m.group(2)); break
except FileNotFoundError:
    pass
if best is not None:
    with open(f"{acc}_dock.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["structure", "smina_affinity_kcal"]); w.writerow([acc, best])
PY
    """
}
