/*
========================================================================================
    MODULE : COMPRESS_FASTQ — Compression des FASTQ conservés en fin de chaîne
========================================================================================
    Les FASTQ circulent NON compressés entre les étapes du pipeline : compresser
    puis décompresser à chaque étape (gzip, 1 seul cœur) était l'étape la plus
    lente. Seuls les fichiers conservés dans results/ sont compressés, ici, avec
    pigz (gzip multi-cœurs, 4 cœurs), dans une tâche séparée qui tourne en
    parallèle d'EMU sans le retarder.

    Utilisé pour les reads sans l'hôte (results/host_removed/).

    Sorties :
    - {nom du FASTQ}.gz → copie compressée, publiée
========================================================================================
*/

process COMPRESS_FASTQ {

    tag "${sample_id}"
    label 'process_compress'

    publishDir "${params.resultsdir}/host_removed", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${fastq}.gz"), emit: fastq_gz

    // LEÇON : -c écrit sur la sortie standard : le FASTQ d'entrée (un lien vers
    // le fichier de l'étape précédente) n'est ni modifié ni supprimé.
    script:
    """
    pigz -p ${task.cpus} -c "${fastq}" > "${fastq}.gz"
    """
}
