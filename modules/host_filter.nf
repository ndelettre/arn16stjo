/*
========================================================================================
    MODULE : HOST_FILTER — Alignement sur l'hôte, comptage et retrait des reads
========================================================================================
    Pour un échantillon :
      1. minimap2 aligne ses reads filtrés sur l'index du génome de l'hôte
      2. un read est compté comme "hôte" si son alignement couvre au moins
         50 % de sa longueur (on évite ainsi de compter un read bactérien qui
         n'aurait qu'un court alignement fortuit sur le génome humain)
      3. les reads hôte sont retirés : le FASTQ restant part dans EMU

    Chaque échantillon est aligné séparément : il enchaîne vers EMU dès que son
    alignement est fini, sans attendre les autres.

    Sorties :
    - {sample_id}_nohost.fastq    → reads sans l'hôte, non compressés (entrée du
                                    sous-échantillonnage ; la copie publiée est
                                    compressée à part par COMPRESS_FASTQ)
    - {sample_id}_host.tsv        → sample | reads_filtered | host_reads | host_pct | reads_nohost
========================================================================================
*/

process HOST_FILTER {

    tag "${sample_id}"
    label 'process_host_map'

    input:
    tuple val(sample_id), path(fastq)
    path host_index     // index minimap2 (.mmi), commun à tous les échantillons

    output:
    tuple val(sample_id), path("${sample_id}_nohost.fastq"),    emit: reads
    tuple val(sample_id), path("${sample_id}_host.tsv"),        emit: stats

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : colonnes PAF utilisées : $1 = nom du read, $2 = longueur du read,
    // $3/$4 = début/fin de la partie alignée du read.
    //
    // LEÇON : la liste des reads hôte est chargée dans un tableau awk (h[])
    // dans le bloc BEGIN avec getline. On n'utilise pas l'astuce classique
    // "NR == FNR" car elle se trompe quand le premier fichier est vide
    // (échantillon sans aucun read hôte).
    // ─────────────────────────────────────────────────────────────────────────
    script:
    """
    set -euo pipefail

    minimap2 -x map-ont --secondary=no -t ${task.cpus} "${host_index}" "${fastq}" > host.paf

    awk -F'\\t' '(\$4 - \$3) >= 0.5 * \$2 { print \$1 }' host.paf | sort -u > host_ids.txt
    rm -f host.paf

    awk -v ids=host_ids.txt -v sample="${sample_id}" -v stats="${sample_id}_host.tsv" '
        BEGIN { while ((getline id < ids) > 0) h[id] = 1 }
        NR % 4 == 1 { total++; keep = !(substr(\$1, 2) in h); if (!keep) host++ }
        keep { print }
        END {
            pct = (total > 0) ? 100 * host / total : 0
            print "sample\\treads_filtered\\thost_reads\\thost_pct\\treads_nohost" > stats
            printf "%s\\t%d\\t%d\\t%.2f\\t%d\\n", sample, total, host, pct, total - host > stats
        }
      ' "${fastq}" > "${sample_id}_nohost.fastq"

    cat "${sample_id}_host.tsv" >&2
    """
}
