/*
========================================================================================
    MODULE : EMU_ABUNDANCE — Classification taxonomique des reads 16S
========================================================================================
    EMU aligne les reads sur la base 16S (minimap2) puis estime l'abondance
    relative de chaque espèce par un algorithme EM (espérance-maximisation).

    Options utilisées :
    --keep-counts            → ajoute la colonne "estimated counts" (nombre de
                               reads estimé par espèce, utile pour les seuils
                               du rapport et le contrôle interne)
    --min-abundance          → seuil interne d'EMU (défaut 0.0001 = 0,01 %)
    --keep-read-assignments  → (optionnel) probabilité d'attribution de chaque read
    --output-unclassified    → (optionnel) FASTQ des reads non classés

    Sorties :
    - {sample_id}_rel-abundance.tsv                  → abondances + lignée + reads estimés
    - {sample_id}_read-assignment-distributions.tsv  → (optionnel)
    - {sample_id}_*.fastq.gz                         → (optionnel) reads non classés
    - {sample_id}_emu.log                            → journal d'EMU
========================================================================================
*/

process EMU_ABUNDANCE {

    tag "${sample_id}"
    label 'process_emu'

    publishDir "${params.resultsdir}/emu/per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path emu_db

    output:
    tuple val(sample_id), path("${sample_id}_rel-abundance.tsv"), emit: abundance
    path "${sample_id}_read-assignment-distributions.tsv",       emit: assignments, optional: true
    path "${sample_id}_*.fastq.gz",                              emit: unclassified, optional: true
    path "${sample_id}_emu.log",                                 emit: log

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : EMU s'arrête en erreur ("0 reads mapped") si aucun read ne
    // s'aligne sur la base 16S (témoin vide, barcode sans reads...). Pour ne
    // pas bloquer tout le run, on écrit dans ce cas une table au même format
    // qu'EMU, où tous les reads sont "unmapped" (script bin/emu_empty_table.py).
    // Toute autre erreur d'EMU arrête bien le process.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    def opt_assign       = params.emu_keep_read_assignments ? '--keep-read-assignments' : ''
    def opt_unclassified = params.emu_output_unclassified   ? '--output-unclassified'   : ''
    """
    set -euo pipefail

    if [[ "${reads}" == *.gz ]]; then
        n_reads=\$(gzip -cd "${reads}" | awk 'END { print int(NR / 4) }')
    else
        n_reads=\$(awk 'END { print int(NR / 4) }' "${reads}")
    fi

    if [ "\$n_reads" -eq 0 ]; then
        echo "Aucun read à classer pour ${sample_id}" > "${sample_id}_emu.log"
        emu_empty_table.py --taxonomy "${emu_db}/taxonomy.tsv" --unmapped 0 \\
            --output "${sample_id}_rel-abundance.tsv"
    elif emu abundance \\
            --type map-ont \\
            --db "${emu_db}" \\
            --min-abundance ${params.emu_min_abundance} \\
            --keep-counts \\
            ${opt_assign} \\
            ${opt_unclassified} \\
            --threads ${task.cpus} \\
            --output-dir . \\
            --output-basename "${sample_id}" \\
            "${reads}" > "${sample_id}_emu.log" 2>&1; then
        echo "EMU terminé pour ${sample_id}" >&2
    elif grep -q "0 reads mapped" "${sample_id}_emu.log"; then
        echo "Aucun read de ${sample_id} ne s'aligne sur la base 16S (\$n_reads reads non classés)" | tee -a "${sample_id}_emu.log" >&2
        emu_empty_table.py --taxonomy "${emu_db}/taxonomy.tsv" --unmapped "\$n_reads" \\
            --output "${sample_id}_rel-abundance.tsv"
    else
        cat "${sample_id}_emu.log" >&2
        exit 1
    fi

    # Les FASTQ de reads non classés (option) sont compressés. On ne vise que
    # les fichiers écrits par EMU : l'entrée {sample_id}_emu_input.fastq (un
    # lien vers l'étape précédente) a le même préfixe et ne doit pas être touchée.
    for fq in "${sample_id}"_unmapped.fastq "${sample_id}"_unclassified_mapped*.fastq "${sample_id}"_filtered_mapped.fastq; do
        if [ -f "\$fq" ] && [ ! -L "\$fq" ]; then gzip "\$fq"; fi
    done

    # Table complémentaire produite par EMU seulement si --min-abundance dépasse
    # son seuil adaptatif : on ne garde que la table principale.
    rm -f "${sample_id}"_rel-abundance-threshold-*.tsv
    """
}
