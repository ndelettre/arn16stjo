# ARN16StJo

Pipeline Nextflow d'analyse du **16S full-length** séquencé sur **Oxford Nanopore** (R10.4.1), avec classification taxonomique par [EMU](https://github.com/treangenlab/emu) — compatible **EPI2ME Desktop**.

> 🚧 Pipeline en cours de développement. Version actuelle : **0.3.0** (QC, retrait de l'hôte, classification EMU, rapport HTML).

## Description

Le pipeline prend en entrée des reads déjà basecallés et démultiplexés (un FASTQ par échantillon, amorces et adaptateurs retirés), par exemple issus du *Microbial Amplicon Barcoding Kit 24 V14*.

### Workflow

```
FASTQ ─┬─ NANOSTAT (reads bruts)
       └─ CHOPPER (longueur + qualité) ─┬─ NANOSTAT (reads filtrés)
                                        └─ HOST_FILTER (alignement sur l'hôte, reads hôte comptés et retirés)
                                           └─ SUBSAMPLE (≤ 20 000 reads) → EMU_ABUNDANCE ─┬─ EMU_COMBINE (espèce / genre / famille)
                                                             └─ SPIKE_IN (contrôle interne Zymo)
                                                                          ↓
                                              REPORT (rapport HTML unique, tous échantillons)

HOST_INDEX (index du génome hôte, une fois au début) ──→ HOST_FILTER ──→ HOST_INDEX_CLEANUP (suppression de l'index)
```

Étapes :

| Étape | Contenu | Statut |
|---|---|---|
| 1 | QC : NanoStat → Chopper → NanoStat | ✅ |
| 2 | Comptage et retrait des reads de l'hôte (minimap2, GRCh38 par défaut) | ✅ |
| 3 | Classification EMU (base rrnDB 5.10 + NCBI mars 2026 par défaut), contrôle interne Zymo | ✅ |
| 4 | Rapport HTML unique (QC, hôte, contrôle interne, abondances) | ✅ |
| 5 | Mode recherche optionnel (diversité alpha/bêta, export phyloseq) | ⏳ |

## Prérequis

- [Nextflow](https://www.nextflow.io/) ≥ 22.10.0 (EPI2ME Desktop fournit le sien)
- [Docker](https://www.docker.com/)
- Génome de l'hôte (GRCh38 par défaut, ~0,9 Go compressé)
- Base EMU (~180 Mo)

### Téléchargement des bases de données

```bash
# Génome humain GRCh38 (version "no-alt analysis set", recommandée pour l'alignement)
sudo mkdir -p /data/databases/host && sudo chown $USER /data/databases/host
wget -O /data/databases/host/GRCh38.fa.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Base EMU : rrnDB v5.10 + NCBI mars 2026 (fournie dans le dépôt GitHub d'EMU)
sudo mkdir -p /data/databases/emu_db && sudo chown $USER /data/databases/emu_db
cd /data/databases/emu_db
wget https://github.com/treangenlab/emu/raw/master/emu_database/emu_db_mar2026.zip
unzip emu_db_mar2026.zip && rm emu_db_mar2026.zip
# → /data/databases/emu_db/emu_db_mar2026/{species_taxid.fasta, taxonomy.tsv}
```

Le génome de l'hôte peut être fourni :
- en **FASTA** : il est indexé une fois au début de chaque run (quelques minutes, ~10-15 Go de RAM). L'index (~7 Go) est écrit dans le dossier de travail de Nextflow, puis **supprimé dès que tous les échantillons sont alignés** ;
- déjà indexé en **`.mmi`** (`minimap2 -x map-ont -d GRCh38.mmi GRCh38.fa.gz`) : il est utilisé tel quel et n'est jamais supprimé.

## Installation

### Via EPI2ME Desktop

1. Ouvrir EPI2ME Desktop
2. **Workflows** → **Import workflow** → `https://github.com/ndelettre/arn16stjo`
3. Ou **Import Local Workflow** → sélectionner ce dossier

> ⚠️ Le nom du dépôt et `manifest.name` doivent rester **en minuscules** : EPI2ME Desktop
> cherche le workflow dans un dossier au nom converti en minuscules, alors que `nextflow pull`
> garde la casse de l'URL. Avec des majuscules, l'import échoue (« Error refreshing workflow config »).

### En ligne de commande

```bash
nextflow run ndelettre/arn16stjo \
    --fastq_dir /chemin/vers/fastq \
    --spike_in spike_in_I
```

## Données d'entrée

Un dossier contenant **un fichier FASTQ par échantillon** (`.fastq`, `.fq`, `.fastq.gz` ou `.fq.gz`).
Le nom de l'échantillon est le nom du fichier sans l'extension :

| Fichier | Échantillon |
|---|---|
| `patient01.fastq.gz` | `patient01` |
| `temoin_neg.fq` | `temoin_neg` |
| `run2.barcode05.fastq` | `run2.barcode05` |

Les témoins (négatif, positif) sont traités comme des échantillons.
Deux fichiers donnant le même nom d'échantillon (ex. `S1.fastq` et `S1.fastq.gz`) arrêtent le pipeline avec un message d'erreur.

## Paramètres

| Paramètre | Défaut | Description |
|---|---|---|
| `--fastq_dir` | *(obligatoire)* | Dossier contenant les FASTQ |
| `--min_length` | `1200` | Longueur minimale des reads (bp) |
| `--max_length` | `1800` | Longueur maximale des reads (bp) |
| `--min_quality` | `10` | Q-score moyen minimal des reads |
| `--host_reference` | `/data/databases/host/GRCh38.fa.gz` | Génome de l'hôte (FASTA ou `.mmi`) ; vide ou `false` = étape ignorée |
| `--emu_db` | `/data/databases/emu_db/emu_db_mar2026` | Dossier de la base EMU |
| `--emu_max_reads` | `20000` | Reads analysés par EMU au maximum, tirés au hasard (0 = tous) |
| `--emu_min_abundance` | `0.0001` | Seuil interne d'EMU (0,01 %) |
| `--emu_keep_read_assignments` | `false` | Table d'attribution read → espèce |
| `--emu_output_unclassified` | `false` | FASTQ des reads non classés |
| `--spike_in` | `aucun` | Contrôle interne : `aucun`, `spike_in_I` (D6320) ou `spike_in_II` (D6321) |
| `--outdir` | `output` | Dossier des rapports |
| `--resultsdir` | `results` | Dossier des résultats |

### Contrôle interne ZymoBIOMICS

| Valeur | Contrôle | Souches attendues |
|---|---|---|
| `spike_in_I` | Spike-in Control I (D6320, High Microbial Load) | *Imtechella halotolerans*, *Allobacillus halotolerans* |
| `spike_in_II` | Spike-in Control II (D6321, Low Microbial Load) | *Truepera radiovictrix*, *I. halotolerans*, *A. halotolerans* |

Ces trois espèces sont présentes dans la base EMU de mars 2026. Le pipeline vérifie au démarrage qu'elles sont dans la base utilisée et affiche un avertissement sinon.

## Sorties

```
results/
├── nanostat/
│   ├── {échantillon}_raw_nanostat.tsv          # QC avant filtrage
│   └── {échantillon}_filtered_nanostat.tsv     # QC après filtrage
├── chopper/
│   └── {échantillon}_chopper.log               # "Kept X reads out of Y reads"
├── host/
│   └── host_summary.tsv                        # reads hôte par échantillon (nombre et %)
├── host_removed/
│   └── {échantillon}_nohost.fastq.gz           # reads filtrés, sans l'hôte (entrée d'EMU)
└── emu/
    ├── emu-combined-{species,genus,family}.tsv         # abondances relatives, tous échantillons
    ├── emu-combined-{species,genus,family}-counts.tsv  # nombres de reads estimés
    ├── emu_reads_summary.tsv                   # reads sans l'hôte / reads analysés par EMU
    ├── spike_in_counts.tsv                     # reads du contrôle interne par échantillon
    └── per_sample/
        ├── {échantillon}_rel-abundance.tsv     # sortie EMU (abondance, lignée, reads estimés)
        ├── {échantillon}_emu.log
        ├── {échantillon}_read-assignment-distributions.tsv   # option
        └── {échantillon}_{unmapped,unclassified_mapped,filtered_mapped}.fastq.gz  # option
output/
├── rapport_16S.html                            # rapport du run (voir ci-dessous)
├── pipeline_report.html                        # ressources utilisées par tâche
└── pipeline_timeline.html                      # chronologie d'exécution
```

Points à connaître :
- **Reads hôte** : chaque échantillon est aligné séparément sur l'index du génome (10 cœurs par alignement, 2 en parallèle, chacun chargeant ~7-8 Go d'index en RAM). Un read est compté comme hôte si son alignement couvre au moins 50 % de sa longueur. L'index construit pendant le run étant supprimé à la fin, une relance avec `-resume` le reconstruit et refait les alignements sur l'hôte.
- **Lignes non classées** : dans les tables EMU, `unmapped` = reads sans alignement sur la base 16S ; `mapped_filtered` et `mapped_unclassified` = reads alignés mais non attribués. Dans les tables `emu-combined-*`, ces reads sont regroupés sur la ligne sans nom de taxon.
- **Échantillon sans read** (barcode vide, tout filtré, tout hôte, aucun read 16S) : le run continue, l'échantillon apparaît avec 0 read classé.
- **Sous-échantillonnage** : EMU devient très lent au-delà de quelques dizaines de milliers de reads (la base 2026 est très redondante : ~30 reads/s avec 4 cœurs, plus de 12 h pour 1,4 million de reads). Au plus `--emu_max_reads` reads (20 000 par défaut) sont tirés au hasard par échantillon, avec une graine fixe (résultats reproductibles). Avec 20 000 reads, un taxon à 1 % représente ~200 reads. Les nombres de reads d'EMU et du contrôle interne portent sur ces reads analysés ; `emu_reads_summary.tsv` donne, pour chaque échantillon, le nombre de reads avant et après sous-échantillonnage.
- **Seuils d'affichage** (≥ 20 reads et ≥ 1 %) : appliqués dans le rapport HTML et modifiables dans la page ; les tables ci-dessus contiennent tous les taxons.

### Rapport HTML (`output/rapport_16S.html`)

Un seul fichier, au style du rapport phyloreport (logo HPSJ), à ouvrir dans un navigateur (Plotly est chargé depuis Internet) :

- **Paramètres du run** et chiffres clés (échantillons, reads bruts, après Chopper, analysés par EMU) ;
- **Résumé par échantillon** : reads bruts, après Chopper, reads hôte (nombre et %), sans l'hôte, analysés et classés par EMU ;
- **Contrôle qualité** (NanoStat) avant / après Chopper : graphique des reads (échelle log) et tableau (reads, Mb, longueur médiane, Q médiane, N50) ;
- **Contrôle interne** : reads estimés et abondance (% des reads classés) de chaque souche, par échantillon ;
- **Abondances relatives** (EMU) :
  - barres empilées par échantillon, au rang **espèce, genre ou famille** ;
  - un taxon est affiché s'il atteint **à la fois** le nombre de reads minimum (20 par défaut) **et** l'abondance minimale (1 % par défaut) ; ces deux seuils se modifient dans la page ; le reste est regroupé dans « Sous les seuils » ;
  - le contrôle interne est retiré et les abondances recalculées sur 100 % des reads classés ;
  - les 8 taxons les plus abondants ont chacun leur couleur (fixe, quels que soient les seuils), les autres sont regroupés dans « Autres taxons » (détail au survol) ;
  - tableau détaillé (échantillon, taxon, reads, %) téléchargeable en **CSV** (séparateur « ; », compatible Excel).

## Images Docker

Toutes les images sont figées sur une version exacte pour garantir la reproductibilité.

| Outil | Image |
|---|---|
| NanoStat | `quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0` |
| Chopper | `quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0` |
| minimap2 (hôte) | `staphb/minimap2:2.31` |
| EMU | `quay.io/biocontainers/emu:3.6.2--hdfd78af_0` |
| Rapport (Python) | `nicolasdelettre/phyloreport:1.0.1` |
