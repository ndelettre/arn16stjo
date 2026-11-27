/*
========================================================================================
    MODULE : SUBSAMPLE — Sous-échantillonnage des reads avant EMU
========================================================================================
    EMU devient très lent au-delà de quelques dizaines de milliers de reads
    (sur la base 2026 : ~30 reads/s avec 4 cœurs, soit plus de 12 h pour un
    échantillon de 1,4 million de reads). On tire donc au hasard au plus
    --emu_max_reads reads par échantillon (20 000 par défaut ; 0 = tous).

    Le tirage utilise une graine fixe (bin/subsample_fastq.py) : relancer le
    pipeline sur les mêmes données donne exactement les mêmes reads.

    Avec 20 000 reads, un taxon à 1 % représente ~200 reads : les abondances
    relatives restent fiables. Les comptes de reads d'EMU (et du contrôle
    interne) portent sur les reads analysés, pas sur la totalité.

    Sorties :
    - {sample_id}_emu_input.fastq    → reads analysés par EMU (non compressé)
    - {sample_id}_subsample.tsv      → sample | reads_nohost | reads_emu | subsampled
========================================================================================
*/

process SUBSAMPLE {

    tag "${sample_id}"
    label 'process_mono'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_emu_input.fastq"),    emit: reads
    tuple val(sample_id), path("${sample_id}_subsample.tsv"),      emit: stats

    script:
    """
    subsample_fastq.py \\
        --input "${fastq}" \\
        --max-reads ${params.emu_max_reads} \\
        --seed 42 \\
        --sample "${sample_id}" \\
        --output "${sample_id}_emu_input.fastq" \\
        --stats "${sample_id}_subsample.tsv"

    cat "${sample_id}_subsample.tsv" >&2
    """
}
