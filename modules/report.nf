/*
========================================================================================
    MODULE : REPORT — Rapport HTML unique du run
========================================================================================
    Rassemble les résultats de tous les échantillons dans un seul fichier HTML
    (script bin/report.py, style phyloreport, logo HPSJ) :
      - résumé par échantillon (reads bruts, après Chopper, hôte, analysés par EMU)
      - contrôle qualité NanoStat avant / après Chopper
      - contrôle interne ZymoBIOMICS (reads et abondance par souche)
      - abondances EMU par espèce / genre / famille : seuils (reads ET %)
        modifiables dans la page, tableau détaillé exportable en CSV

    Sortie :
    - rapport_16S.html (dans --outdir)
========================================================================================
*/

// Entoure une valeur d'apostrophes pour le shell (chemins avec espaces...)
// en échappant les apostrophes qu'elle pourrait contenir.
def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process REPORT {

    label 'process_mono'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path result_files       // TSV NanoStat, hôte, sous-échantillonnage, EMU (tous échantillons)
    path logo
    val spike_name          // nom du contrôle interne ('' si aucun)
    val spike_species       // espèces du contrôle, séparées par des virgules ('' si aucun)

    output:
    path "rapport_16S.html", emit: html

    // ─────────────────────────────────────────────────────────────────────────
    // LEÇON : une fonction (def shellQuote) déclarée dans le module peut être
    // utilisée dans le script du process.
    // ─────────────────────────────────────────────────────────────────────────
    script:
    def run_params = [
        "Dossier FASTQ=${params.fastq_dir}",
        "Filtrage Chopper=${params.min_length}–${params.max_length} bp, Q ≥ ${params.min_quality}",
        "Génome hôte=${params.host_reference ?: 'aucun'}",
        "Base EMU=${params.emu_db}",
        "Reads analysés par EMU=${params.emu_max_reads ? "au plus ${params.emu_max_reads} par échantillon" : 'tous'}",
        "Seuil interne EMU=${params.emu_min_abundance}",
        "Contrôle interne=${spike_name ?: 'aucun'}",
    ].collect { item -> "--param ${shellQuote(item)}" }.join(' ')
    def tools = [
        "NanoStat=1.6.0", "Chopper=0.9.0", "minimap2=2.31", "EMU=3.6.2", "Nextflow=${nextflow.version}",
    ].collect { item -> "--tool ${shellQuote(item)}" }.join(' ')
    """
    report.py \\
        --input-dir . \\
        --output rapport_16S.html \\
        --logo "${logo}" \\
        --version ${shellQuote(workflow.manifest.version)} \\
        --spike-name ${shellQuote(spike_name)} \\
        --spike-species ${shellQuote(spike_species)} \\
        --max-reads ${params.emu_max_reads} \\
        --min-reads 20 \\
        --min-pct 1 \\
        ${run_params} \\
        ${tools}
    """
}
