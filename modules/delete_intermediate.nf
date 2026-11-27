/*
========================================================================================
    MODULE : DELETE_INTERMEDIATE — Suppression d'un FASTQ intermédiaire
========================================================================================
    Les FASTQ circulent non compressés entre les étapes : pour un run de 24
    échantillons profonds, chaque série de fichiers intermédiaires pèse
    ~80 Go dans le dossier work/. Ce module supprime un FASTQ intermédiaire
    dès que toutes les étapes qui le lisent sont terminées (main.nf attend
    leurs résultats avec .join() avant de l'appeler).

    Utilisé (alias dans main.nf) pour :
    - CLEANUP_CHOPPER : sortie de Chopper, une fois lue par NanoStat et par
                        le retrait de l'hôte (ou par EMU s'il n'y a pas d'hôte)
    - CLEANUP_NOHOST  : reads sans l'hôte, une fois compressés (COMPRESS_FASTQ),
                        comptés/tirés (SUBSAMPLE) et classés (EMU)

    Conséquence pour "-resume" : les fichiers supprimés sont recalculés si le
    run est relancé, et les étapes suivantes sont refaites.
========================================================================================
*/

process DELETE_INTERMEDIATE {

    tag "${sample_id}"

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : comme HOST_INDEX_CLEANUP, un process "exec:" (code Groovy exécuté
    // par Nextflow, sans conteneur). Le fichier est reçu en "val" pour avoir
    // son chemin réel dans le dossier de travail de l'étape qui l'a produit.
    // ─────────────────────────────────────────────────────────────────────────
    input:
    tuple val(sample_id), val(fastq)

    exec:
    def fastq_file = fastq.toRealPath()
    def size_gb = String.format('%.2f', fastq_file.size() / 1e9)
    fastq_file.delete()
    log.info "${sample_id} : fichier intermédiaire supprimé (${size_gb} Go) : ${fastq_file.name}"
}
