#!/usr/bin/env nextflow
/*
========================================================================================
    ARN16StJo — PIPELINE 16S FULL-LENGTH NANOPORE (EMU)  v0.3.0
    Compatible EPI2ME | Nextflow DSL2
========================================================================================
    WORKFLOW :
       NANOSTAT (reads bruts)
       CHOPPER → NANOSTAT (reads filtrés)
               → HOST_FILTER (alignement sur l'hôte, comptage et retrait des reads)
                 (index du génome : HOST_INDEX au début, HOST_INDEX_CLEANUP à la fin)
               → SUBSAMPLE (au plus --emu_max_reads reads par échantillon)
               → EMU_ABUNDANCE → EMU_COMBINE (tables espèce / genre / famille)
                               → SPIKE_IN (reads du contrôle interne Zymo)
       → REPORT (rapport HTML unique : QC, hôte, contrôle interne, abondances)

    Étape à venir : mode recherche (diversité).
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// IMPORTS DES MODULES
// LEÇON : un même process peut être importé deux fois sous deux noms
// différents avec "as". Chaque alias est un process indépendant : on peut
// donc lancer NanoStat avant ET après Chopper sans dupliquer le module.
// ─────────────────────────────────────────────────────────────────────────────
include { NANOSTAT as NANOSTAT_RAW      } from './modules/nanostat.nf' // QC reads bruts
include { CHOPPER                       } from './modules/chopper.nf'  // Filtrage longueur/qualité
include { NANOSTAT as NANOSTAT_FILTERED } from './modules/nanostat.nf' // QC reads filtrés
include { HOST_INDEX                    } from './modules/host_index.nf'         // Index du génome hôte
include { HOST_FILTER                   } from './modules/host_filter.nf'        // Alignement/retrait hôte
include { HOST_INDEX_CLEANUP            } from './modules/host_index_cleanup.nf' // Suppression de l'index
include { SUBSAMPLE                     } from './modules/subsample.nf'     // Sous-échantillonnage
include { EMU_ABUNDANCE                 } from './modules/emu_abundance.nf' // Classification EMU
include { EMU_COMBINE                   } from './modules/emu_combine.nf'   // Tables multi-échantillons
include { SPIKE_IN                      } from './modules/spike_in.nf'      // Contrôle interne Zymo
include { REPORT                        } from './modules/report.nf'        // Rapport HTML

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW PRINCIPAL
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // LEÇON : la bannière et les validations sont DANS le bloc workflow :
    // la syntaxe stricte de Nextflow (par défaut dans les futures versions)
    // interdit les instructions en dehors d'un workflow ou d'un process.

    // ── Génome hôte (optionnel) ───────────────────────────────────────────────
    // LEÇON : un paramètre laissé vide arrive comme "" (EPI2ME) ou comme true
    // (ligne de commande "--host_reference" sans valeur). On ne garde le chemin
    // que s'il s'agit d'un texte non vide ; sinon l'étape hôte est ignorée.
    def host_reference = (params.host_reference instanceof CharSequence && params.host_reference.toString().trim())
        ? params.host_reference.toString().trim()
        : null

    // ── Bannière de démarrage ─────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║     ARN16StJo - 16S FULL-LENGTH NANOPORE (EMU)  v0.3.0   ║
    ╚══════════════════════════════════════════════════════════╝
      Dossier FASTQ   : ${params.fastq_dir}
      Dossier sortie  : ${params.outdir}
      Longueur reads  : ${params.min_length} - ${params.max_length} bp
      Qualité min.    : Q${params.min_quality}
      Génome hôte     : ${host_reference ?: 'aucun (étape ignorée)'}
      Base EMU        : ${params.emu_db}
      Reads max. EMU  : ${params.emu_max_reads ?: 'tous'}
      Contrôle interne: ${params.spike_in}
    ──────────────────────────────────────────────────────────
    """.stripIndent()

    // ── Validation des paramètres ─────────────────────────────────────────────
    if (!params.fastq_dir) {
        error "ERREUR : --fastq_dir est obligatoire.\nUsage : nextflow run main.nf --fastq_dir /chemin/vers/fastq"
    }
    if (params.min_length >= params.max_length) {
        error "ERREUR : --min_length (${params.min_length}) doit être inférieur à --max_length (${params.max_length})."
    }
    if (host_reference && !(host_reference =~ /\.((fa|fasta|fna)(\.gz)?|mmi)$/)) {
        error "ERREUR : --host_reference doit être un FASTA (.fa, .fasta, .fna, éventuellement .gz) ou un index minimap2 (.mmi)."
    }
    if (host_reference && !file(host_reference).exists()) {
        error "ERREUR : génome hôte introuvable : ${host_reference}\n(voir le README pour le télécharger, ou laisser --host_reference vide pour ignorer l'étape)."
    }
    if (!(params.emu_max_reads instanceof Integer) || params.emu_max_reads < 0) {
        error "ERREUR : --emu_max_reads doit être un entier ≥ 0 (0 = pas de sous-échantillonnage ; reçu : '${params.emu_max_reads}')."
    }
    ['species_taxid.fasta', 'taxonomy.tsv'].each { db_file ->
        if (!file("${params.emu_db}/${db_file}").exists()) {
            error "ERREUR : ${db_file} introuvable dans la base EMU --emu_db (${params.emu_db}).\n(voir le README pour télécharger la base)."
        }
    }

    // ── Contrôle interne ZymoBIOMICS ──────────────────────────────────────────
    // Souches attendues pour chaque contrôle (fiches Zymo D6320 et D6321).
    def spike_in_controls = [
        spike_in_I : [name   : 'Spike-in Control I (D6320, High Microbial Load)',
                      species: ['Imtechella halotolerans', 'Allobacillus halotolerans']],
        spike_in_II: [name   : 'Spike-in Control II (D6321, Low Microbial Load)',
                      species: ['Truepera radiovictrix', 'Imtechella halotolerans', 'Allobacillus halotolerans']]
    ]
    if (params.spike_in != 'aucun' && !spike_in_controls.containsKey(params.spike_in)) {
        error "ERREUR : --spike_in doit valoir 'aucun', 'spike_in_I' ou 'spike_in_II' (reçu : '${params.spike_in}')."
    }
    def spike_in = spike_in_controls[params.spike_in]

    // Vérifie que les souches du contrôle existent dans la base EMU : sinon,
    // leurs reads seraient attribués à une espèce proche et fausseraient tout.
    if (spike_in) {
        def db_species = file("${params.emu_db}/taxonomy.tsv").readLines().collect { line -> line.tokenize('\t')[1] } as Set
        spike_in.species.findAll { sp -> !db_species.contains(sp) }.each { sp ->
            log.warn "Contrôle interne : '${sp}' est absente de la base EMU ${params.emu_db}"
        }
    }

    // ── Création du channel d'entrée ──────────────────────────────────────────
    // Un fichier FASTQ = un échantillon. Le nom de l'échantillon est le nom du
    // fichier sans l'extension .fastq/.fq(.gz) : "patient.01.fastq.gz" → "patient.01".
    // LEÇON : on n'utilise pas getSimpleName() car il coupe au PREMIER point
    // ("patient.01.fastq.gz" deviendrait "patient").
    //
    // LEÇON : groupTuple() regroupe les fichiers par nom d'échantillon. Si deux
    // fichiers donnent le même nom (ex : S1.fastq et S1.fastq.gz), on arrête
    // le pipeline au lieu d'écraser silencieusement un des deux.
    ch_fastq = channel
        .fromPath(
            "${params.fastq_dir}/*.{fastq,fastq.gz,fq,fq.gz}",
            checkIfExists: true
        )
        .map { file ->
            def sample_id = file.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample_id, file)
        }
        .groupTuple()
        .map { sample_id, files ->
            if (files.size() > 1) {
                error "ERREUR : plusieurs fichiers pour l'échantillon '${sample_id}' : ${files*.name.join(', ')}"
            }
            tuple(sample_id, files[0])
        }

    // ── Contrôle qualité avant filtrage ───────────────────────────────────────
    // LEÇON : le deuxième argument ("raw") est une simple valeur passée au
    // process : elle sert à nommer les fichiers et le dossier de sortie.
    NANOSTAT_RAW(ch_fastq, "raw")

    // ── Filtrage longueur / qualité ───────────────────────────────────────────
    // Chopper garde les reads entre min_length et max_length avec un Q-score
    // moyen ≥ min_quality (16S full-length ≈ 1 500 bp).
    CHOPPER(ch_fastq)

    // Alerte dans le log EPI2ME pour chaque échantillon sans read après filtrage
    // (le run continue : l'échantillon apparaîtra avec 0 read dans le rapport).
    CHOPPER.out.log
        .filter { _sample_id, chopper_log -> chopper_log.text =~ /^Kept 0 reads/ }
        .subscribe { sample_id, chopper_log ->
            log.warn "${sample_id} : aucun read conservé par Chopper (${chopper_log.text.trim()})"
        }

    // ── Contrôle qualité après filtrage ───────────────────────────────────────
    NANOSTAT_FILTERED(CHOPPER.out.reads, "filtered")

    // ── Reads de l'hôte : comptage et retrait ─────────────────────────────────
    // Le génome de l'hôte est indexé UNE fois au début du run (HOST_INDEX,
    // en parallèle de Chopper), puis chaque échantillon est aligné séparément
    // sur cet index (HOST_FILTER). L'index est supprimé dès que tous les
    // échantillons sont alignés (HOST_INDEX_CLEANUP). Un index .mmi fourni
    // par l'utilisateur est utilisé tel quel et n'est jamais supprimé.
    //
    // LEÇON : .first() transforme le channel en "value channel" : l'index
    // est réutilisé par HOST_FILTER pour chacun des échantillons au lieu
    // d'être consommé par le premier.
    if (host_reference) {
        ch_host_reference = channel.fromPath(host_reference, checkIfExists: true).first()

        if (host_reference.endsWith('.mmi')) {
            ch_host_index = ch_host_reference
        } else {
            HOST_INDEX(ch_host_reference)
            ch_host_index = HOST_INDEX.out.index.first()
        }

        HOST_FILTER(CHOPPER.out.reads, ch_host_index)
        ch_reads_emu   = HOST_FILTER.out.reads
        ch_host_stats  = HOST_FILTER.out.stats.map { _sample_id, stats -> stats }

        if (!host_reference.endsWith('.mmi')) {
            HOST_INDEX_CLEANUP(ch_host_index, HOST_FILTER.out.stats.collect())
        }

        // Table récapitulative : une ligne par échantillon
        HOST_FILTER.out.stats
            .map { _sample_id, stats -> stats }
            .collectFile(name: 'host_summary.tsv', keepHeader: true, skip: 1,
                         sort: { stats -> stats.name }, storeDir: "${params.resultsdir}/host")
    } else {
        log.warn "Aucun génome hôte (--host_reference) : les reads de l'hôte ne sont ni comptés ni retirés."
        ch_reads_emu  = CHOPPER.out.reads
        ch_host_stats = channel.empty()
    }

    // ── Sous-échantillonnage avant EMU ────────────────────────────────────────
    // EMU est très lent au-delà de quelques dizaines de milliers de reads :
    // on en tire au hasard au plus --emu_max_reads par échantillon (0 = tous).
    SUBSAMPLE(ch_reads_emu)

    SUBSAMPLE.out.stats
        .map { _sample_id, stats -> stats }
        .collectFile(name: 'emu_reads_summary.tsv', keepHeader: true, skip: 1,
                     sort: { stats -> stats.name }, storeDir: "${params.resultsdir}/emu")

    // ── Classification taxonomique EMU ────────────────────────────────────────
    ch_emu_db = channel.fromPath(params.emu_db, type: 'dir', checkIfExists: true).first()
    EMU_ABUNDANCE(SUBSAMPLE.out.reads, ch_emu_db)

    ch_abundance_tables = EMU_ABUNDANCE.out.abundance
        .map { _sample_id, table -> table }
        .collect()

    // Tables multi-échantillons (espèce, genre, famille)
    EMU_COMBINE(ch_abundance_tables)

    // ── Contrôle interne ──────────────────────────────────────────────────────
    if (spike_in) {
        SPIKE_IN(ch_abundance_tables, spike_in.name, spike_in.species.join(','))
    }

    // ── Rapport HTML ──────────────────────────────────────────────────────────
    // LEÇON : .mix() fusionne les fichiers de résultats de toutes les étapes,
    // .collect() attend qu'ils soient tous produits : le rapport est généré
    // une seule fois, à la fin, avec tous les échantillons.
    ch_report_files = NANOSTAT_RAW.out.stats.map { _sample_id, stats -> stats }
        .mix(NANOSTAT_FILTERED.out.stats.map { _sample_id, stats -> stats })
        .mix(ch_host_stats)
        .mix(SUBSAMPLE.out.stats.map { _sample_id, stats -> stats })
        .mix(EMU_ABUNDANCE.out.abundance.map { _sample_id, table -> table })
        .collect()

    REPORT(
        ch_report_files,
        file("${projectDir}/assets/logo_hpsj.webp", checkIfExists: true),
        spike_in ? spike_in.name : '',
        spike_in ? spike_in.species.join(',') : ''
    )

    // ── Résumé de fin de pipeline (succès ou échec) ───────────────────────────
    // LEÇON : à l'intérieur du bloc onComplete, les noms `workflow` et `params`
    // ne sont plus accessibles : on en garde des copies locales avant.
    def run        = workflow
    def outdir     = params.outdir
    def resultsdir = params.resultsdir
    run.onComplete {
        log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║                  PIPELINE TERMINÉ !                      ║
    ╚══════════════════════════════════════════════════════════╝
      Statut    : ${run.success ? '✅ Succès' : '❌ Échec'}
      Durée     : ${run.duration}
      Rapports  : ${outdir}/
      Résultats : ${resultsdir}/
    ──────────────────────────────────────────────────────────
    """.stripIndent()
    }
}
