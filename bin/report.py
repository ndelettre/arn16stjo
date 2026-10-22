#!/usr/bin/env python3
"""
report.py
---------
Génère le rapport HTML unique du pipeline ARN16StJo (style phyloreport) :

  - En-tête (logo HPSJ, date, version) et paramètres du run
  - Résumé par échantillon (reads bruts, après Chopper, hôte, analysés par EMU)
  - Contrôle qualité NanoStat avant / après Chopper (tableau + graphique)
  - Contrôle interne ZymoBIOMICS (reads et abondance par souche)
  - Abondances EMU : barres empilées par espèce / genre / famille, seuils
    (reads minimum ET abondance minimale) modifiables dans la page, tableau
    détaillé exportable en CSV. Le contrôle interne est retiré et les
    abondances recalculées sur 100 %.

Le script n'utilise que la bibliothèque standard de Python (3.7+) : les
graphiques sont dessinés dans le navigateur par Plotly.js (chargé depuis
Internet, comme pour phyloreport).

Fichiers lus dans --input-dir :
    {sample}_raw_nanostat.tsv, {sample}_filtered_nanostat.tsv   (NanoStat --tsv)
    {sample}_host.tsv                                           (HOST_FILTER, optionnel)
    {sample}_subsample.tsv                                      (SUBSAMPLE)
    {sample}_rel-abundance.tsv                                  (EMU)

Usage :
    report.py --input-dir . --output rapport_16S.html --logo logo.webp \
        --version 0.3.0 --spike-name "Spike-in Control II (D6321)" \
        --spike-species "Truepera radiovictrix,Imtechella halotolerans,..." \
        --param "Base EMU=/data/databases/emu_db/emu_db_mar2026" --tool "EMU=3.6.2"
"""

import argparse
import base64
import csv
import glob
import html
import json
import os
import re
from datetime import datetime

NOT_CLASSIFIED = {"unmapped", "mapped_filtered", "mapped_unclassified"}
RANKS = [("species", "Espèce"), ("genus", "Genre"), ("family", "Famille")]


# ─────────────────────────────────────────────────────────────────────────────
# Lecture des fichiers
# ─────────────────────────────────────────────────────────────────────────────
def natural_key(text):
    """Tri naturel : barcode2 avant barcode10."""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def files_by_sample(input_dir, suffix):
    found = {}
    for path in glob.glob(os.path.join(input_dir, "*" + suffix)):
        found[os.path.basename(path)[: -len(suffix)]] = path
    return found


def read_nanostat(path):
    """NanoStat --tsv : lignes 'métrique<TAB>valeur'."""
    stats = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "Metrics":
                stats[parts[0]] = parts[1]

    def num(key):
        try:
            return float(stats.get(key, 0))
        except ValueError:
            return 0.0

    return {
        "reads": int(num("number_of_reads")),
        "bases": num("number_of_bases"),
        "median_length": num("median_read_length"),
        "median_qual": num("median_qual"),
        "n50": num("n50"),
    }


def read_single_row_tsv(path):
    with open(path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0] if rows else {}


def read_emu(path, spike_species):
    """
    Table EMU {sample}_rel-abundance.tsv (avec "estimated counts").
    Renvoie les reads estimés par taxon (hors contrôle interne) pour chaque
    rang, et les totaux utiles au résumé.
    """
    by_rank = {rank: {} for rank, _ in RANKS}
    spike = {}
    total = classified = 0.0
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                reads = float(row.get("estimated counts") or 0)
            except ValueError:
                reads = 0.0
            total += reads
            if row.get("tax_id") in NOT_CLASSIFIED:
                continue
            classified += reads
            species = (row.get("species") or "").strip()
            if species in spike_species:
                spike[species] = spike.get(species, 0.0) + reads
                continue
            for rank, label in RANKS:
                name = (row.get(rank) or "").strip() or "Non attribué ({})".format(label.lower())
                by_rank[rank][name] = by_rank[rank].get(name, 0.0) + reads
    return {
        "by_rank": by_rank,
        "spike": spike,
        "total": total,
        "classified": classified,
        "bacteria": classified - sum(spike.values()),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Mise en forme
# ─────────────────────────────────────────────────────────────────────────────
def fmt_int(value):
    if value is None:
        return '<span class="dash">—</span>'
    return "{:,}".format(int(round(value))).replace(",", " ")


def fmt_float(value, digits=1):
    if value is None:
        return '<span class="dash">—</span>'
    return "{:.{d}f}".format(value, d=digits).replace(".", ",")


def fmt_pct(part, whole, digits=1):
    if not whole:
        return '<span class="dash">—</span>'
    return fmt_float(100.0 * part / whole, digits) + " %"


def esc(text):
    return html.escape(str(text))


def section(title, icon_path, body):
    return """
    <div class="section">
        <div class="section-header">
            <svg class="section-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="{icon}"/></svg>
            <h2>{title}</h2>
        </div>
        <div class="section-body">{body}</div>
    </div>""".format(icon=icon_path, title=esc(title), body=body)


ICON_LIST = "M4 6h16M4 12h16M4 18h16"
ICON_CHECK = ("M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04"
              "A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 "
              "0-1.042-.133-2.052-.382-3.016z")
ICON_VIAL = "M9 3h6M10 3v6l-5 9a2 2 0 002 3h10a2 2 0 002-3l-5-9V3"
ICON_GEAR = ("M12 15a3 3 0 100-6 3 3 0 000 6zM19.4 15a1.65 1.65 0 00.33 1.82l.06.06a2 2 0 11-2.83 "
             "2.83l-.06-.06a1.65 1.65 0 00-1.82-.33 1.65 1.65 0 00-1 1.51V21a2 2 0 11-4 0v-.09"
             "A1.65 1.65 0 009 19.4a1.65 1.65 0 00-1.82.33l-.06.06a2 2 0 11-2.83-2.83l.06-.06"
             "A1.65 1.65 0 004.6 15a1.65 1.65 0 00-1.51-1H3a2 2 0 110-4h.09A1.65 1.65 0 004.6 9"
             "a1.65 1.65 0 00-.33-1.82l-.06-.06a2 2 0 112.83-2.83l.06.06A1.65 1.65 0 009 4.6"
             "a1.65 1.65 0 001-1.51V3a2 2 0 114 0v.09a1.65 1.65 0 001 1.51 1.65 1.65 0 001.82-.33"
             "l.06-.06a2 2 0 112.83 2.83l-.06.06A1.65 1.65 0 0019.4 9a1.65 1.65 0 001.51 1H21"
             "a2 2 0 110 4h-.09a1.65 1.65 0 00-1.51 1z")


def table(headers, rows, numeric_from=1):
    head = "".join(
        '<th class="num">{}</th>'.format(h) if i >= numeric_from else "<th>{}</th>".format(h)
        for i, h in enumerate(headers))
    body = "".join(
        "<tr>" + "".join(
            '<td class="num">{}</td>'.format(c) if i >= numeric_from else "<td>{}</td>".format(c)
            for i, c in enumerate(row)) + "</tr>"
        for row in rows)
    return '<div class="table-wrap"><table><thead><tr>{}</tr></thead><tbody>{}</tbody></table></div>'.format(head, body)


# ─────────────────────────────────────────────────────────────────────────────
# Rapport
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", default=".")
    parser.add_argument("--output", required=True)
    parser.add_argument("--logo", help="logo (WEBP) affiché dans l'en-tête")
    parser.add_argument("--version", default="")
    parser.add_argument("--title", default="Rapport d'analyse 16S full-length")
    parser.add_argument("--spike-name", default="")
    parser.add_argument("--spike-species", default="", help="espèces du contrôle, séparées par des virgules")
    parser.add_argument("--max-reads", type=int, default=0, help="--emu_max_reads utilisé")
    parser.add_argument("--min-reads", type=float, default=20, help="seuil d'affichage par défaut (reads)")
    parser.add_argument("--min-pct", type=float, default=1, help="seuil d'affichage par défaut (%%)")
    parser.add_argument("--param", action="append", default=[], help="'Nom=valeur' affiché dans les paramètres")
    parser.add_argument("--tool", action="append", default=[], help="'Outil=version' affiché en pied de page")
    args = parser.parse_args()

    spike_species = [s.strip() for s in args.spike_species.split(",") if s.strip()]

    raw = {s: read_nanostat(p) for s, p in files_by_sample(args.input_dir, "_raw_nanostat.tsv").items()}
    filt = {s: read_nanostat(p) for s, p in files_by_sample(args.input_dir, "_filtered_nanostat.tsv").items()}
    host = {s: read_single_row_tsv(p) for s, p in files_by_sample(args.input_dir, "_host.tsv").items()}
    sub = {s: read_single_row_tsv(p) for s, p in files_by_sample(args.input_dir, "_subsample.tsv").items()}
    emu = {s: read_emu(p, spike_species) for s, p in files_by_sample(args.input_dir, "_rel-abundance.tsv").items()}

    samples = sorted(set(raw) | set(filt) | set(emu), key=natural_key)

    # ── Résumé par échantillon ───────────────────────────────────────────────
    summary_rows = []
    for s in samples:
        r, f, h, sb, e = raw.get(s), filt.get(s), host.get(s), sub.get(s), emu.get(s)
        host_reads = int(h["host_reads"]) if h else None
        reads_filtered = f["reads"] if f else None
        summary_rows.append([
            esc(s),
            fmt_int(r["reads"] if r else None),
            fmt_int(reads_filtered),
            (fmt_int(host_reads) + ' <span class="sub">(' + fmt_pct(host_reads, reads_filtered) + ")</span>")
            if h else '<span class="dash">—</span>',
            fmt_int(int(sb["reads_nohost"]) if sb else None),
            fmt_int(int(sb["reads_emu"]) if sb else None),
            (fmt_int(e["classified"]) + ' <span class="sub">(' + fmt_pct(e["classified"], e["total"]) + ")</span>")
            if e else '<span class="dash">—</span>',
        ])
    summary_html = table(
        ["Échantillon", "Reads bruts", "Après Chopper", "Reads hôte", "Sans l'hôte",
         "Analysés par EMU", "Classés par EMU"], summary_rows)
    if args.max_reads:
        summary_html += ('<p class="note">EMU analyse au plus {} reads par échantillon, tirés au hasard '
                         "(graine fixe). Les abondances et le contrôle interne portent sur ces reads.</p>"
                         ).format(fmt_int(args.max_reads))

    # ── Contrôle qualité ─────────────────────────────────────────────────────
    qc_rows = []
    for s in samples:
        r, f = raw.get(s), filt.get(s)
        row = [esc(s)]
        for d in (r, f):
            row += [fmt_int(d["reads"]) if d else "—",
                    fmt_float(d["bases"] / 1e6, 1) if d else "—",
                    fmt_int(d["median_length"]) if d else "—",
                    fmt_float(d["median_qual"], 1) if d else "—",
                    fmt_int(d["n50"]) if d else "—"]
        qc_rows.append(row)
    qc_head = ('<div class="table-wrap"><table><thead>'
               '<tr><th rowspan="2">Échantillon</th><th class="grp" colspan="5">Avant Chopper</th>'
               '<th class="grp" colspan="5">Après Chopper</th></tr><tr>'
               + "".join('<th class="num">{}</th>'.format(h) for h in
                         ["Reads", "Mb", "Long. médiane", "Q médiane", "N50"] * 2)
               + "</tr></thead><tbody>"
               + "".join("<tr>" + "".join(
                   ('<td class="num">{}</td>' if i else "<td>{}</td>").format(c) for i, c in enumerate(row))
                   + "</tr>" for row in qc_rows)
               + "</tbody></table></div>")
    qc_html = '<div id="qc_chart" class="chart"></div>' + qc_head

    # ── Contrôle interne ─────────────────────────────────────────────────────
    if spike_species:
        spike_rows = []
        for s in samples:
            e = emu.get(s)
            row = [esc(s)]
            for sp in spike_species:
                if e is None:
                    row.append('<span class="dash">—</span>')
                    continue
                reads = e["spike"].get(sp, 0.0)
                row.append(fmt_int(reads) + ' <span class="sub">(' + fmt_pct(reads, e["classified"], 2) + ")</span>")
            spike_rows.append(row)
        spike_html = ('<p class="lead">{}</p>'.format(esc(args.spike_name))
                      + table(["Échantillon"] + ['<i>{}</i>'.format(esc(sp)) for sp in spike_species], spike_rows)
                      + '<p class="note">Reads estimés par EMU (abondance en % des reads classés). '
                        "Le contrôle interne est retiré des graphiques et tableaux d'abondance.</p>")
    else:
        spike_html = '<p class="lead">Aucun contrôle interne déclaré pour ce run.</p>'

    # ── Données des abondances pour le graphique interactif ─────────────────
    abundance = {"samples": samples, "totals": {}, "classified": {}, "ranks": {}, "order": {},
                 "minReads": args.min_reads, "minPct": args.min_pct}
    for s in samples:
        abundance["totals"][s] = round(emu[s]["bacteria"], 3) if s in emu else 0
        abundance["classified"][s] = round(emu[s]["classified"], 3) if s in emu else 0
    for rank, _ in RANKS:
        taxa = {}
        for s in samples:
            if s not in emu:
                continue
            for name, reads in emu[s]["by_rank"][rank].items():
                taxa.setdefault(name, {})[s] = round(reads, 3)
        # Ordre fixe des couleurs : taxons classés par leur abondance maximale
        # dans un échantillon (un pathogène dominant d'un seul échantillon
        # garde sa couleur), puis par abondance cumulée.
        def share(name):
            values = [reads / abundance["totals"][s] for s, reads in taxa[name].items()
                      if abundance["totals"][s] > 0]
            return (max(values) if values else 0, sum(values))
        abundance["ranks"][rank] = taxa
        abundance["order"][rank] = sorted(taxa, key=share, reverse=True)

    qc_data = {"samples": samples,
               "raw": [raw[s]["reads"] if s in raw else None for s in samples],
               "filtered": [filt[s]["reads"] if s in filt else None for s in samples]}

    # ── Paramètres et pied de page ───────────────────────────────────────────
    params_rows = []
    for item in args.param:
        key, _, value = item.partition("=")
        params_rows.append("<tr><th>{}</th><td>{}</td></tr>".format(esc(key), esc(value)))
    params_html = '<table class="kv">{}</table>'.format("".join(params_rows))
    tools_html = "".join('<span class="soft-tag">{} {}</span>'.format(esc(k), esc(v))
                         for k, _, v in (t.partition("=") for t in args.tool))

    logo_html = ""
    if args.logo and os.path.exists(args.logo):
        with open(args.logo, "rb") as handle:
            logo_b64 = base64.b64encode(handle.read()).decode("ascii")
        logo_html = ('<div class="logo"><img src="data:image/webp;base64,{}" '
                     'alt="Hôpitaux Paris Saint-Joseph Marie-Lannelongue"></div><div class="divider"></div>'
                     ).format(logo_b64)

    total_raw = sum(raw[s]["reads"] for s in raw)
    total_filt = sum(filt[s]["reads"] for s in filt)
    total_emu = sum(int(sub[s]["reads_emu"]) for s in sub)
    now = datetime.now()

    def js_json(obj):
        return json.dumps(obj, ensure_ascii=False).replace("</", "<\\/")

    replacements = {
        "%%TITLE%%": esc(args.title),
        "%%LOGO%%": logo_html,
        "%%N_SAMPLES%%": str(len(samples)),
        "%%DATE%%": now.strftime("%d/%m/%Y"),
        "%%VERSION%%": esc("ARN16StJo v" + args.version if args.version else "ARN16StJo"),
        "%%TOTAL_RAW%%": fmt_int(total_raw),
        "%%TOTAL_FILT%%": fmt_int(total_filt),
        "%%PCT_FILT%%": fmt_pct(total_filt, total_raw),
        "%%TOTAL_EMU%%": fmt_int(total_emu),
        "%%PARAMS%%": section("Paramètres du run", ICON_GEAR, params_html),
        "%%SUMMARY%%": section("Résumé par échantillon", ICON_LIST, summary_html),
        "%%QC%%": section("Contrôle qualité des reads (NanoStat)", ICON_CHECK, qc_html),
        "%%SPIKE%%": section("Contrôle interne", ICON_VIAL, spike_html),
        "%%TOOLS%%": tools_html,
        "%%GENERATED%%": now.strftime("%d/%m/%Y à %H:%M"),
        "%%QC_DATA%%": js_json(qc_data),
        "%%ABUNDANCE_DATA%%": js_json(abundance),
    }
    page = TEMPLATE
    for key, value in replacements.items():
        page = page.replace(key, value)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write(page)


# ─────────────────────────────────────────────────────────────────────────────
# Gabarit HTML (style phyloreport). Les graphiques et les seuils sont gérés
# en JavaScript : les seuils se modifient directement dans la page.
# ─────────────────────────────────────────────────────────────────────────────
TEMPLATE = r"""<!DOCTYPE html>
<html lang="fr">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>%%TITLE%%</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f5f6fa; color: #2c3e50; }
    header { background: white; border-bottom: 1px solid #e8eaf0; padding: 1.25rem 2.5rem; display: flex; align-items: center; justify-content: space-between; gap: 1rem; flex-wrap: wrap; }
    .header-left { display: flex; align-items: center; gap: 1.5rem; }
    .logo img { height: 56px; width: auto; }
    .divider { width: 1px; height: 52px; background: #e8eaf0; }
    .header-title h1 { font-size: 18px; font-weight: 600; color: #1a252f; }
    .header-title p { font-size: 13px; color: #7f8c8d; margin-top: 3px; }
    .header-right { text-align: right; }
    .header-right .date { font-size: 14px; font-weight: 500; color: #2c3e50; }
    .header-right .version { font-size: 12px; color: #95a5a6; margin-top: 2px; }
    .container { max-width: 1280px; margin: 0 auto; padding: 1.5rem 2.5rem 2rem; }
    .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; margin-bottom: 1.5rem; }
    .stat-card { background: white; border-radius: 10px; padding: 1rem 1.25rem; border: 1px solid #e8eaf0; }
    .stat-card .value { font-size: 26px; font-weight: 700; color: #2980b9; }
    .stat-card .value small { font-size: 14px; font-weight: 500; color: #7f8c8d; }
    .stat-card .label { font-size: 11px; color: #95a5a6; margin-top: 4px; text-transform: uppercase; letter-spacing: 0.05em; }
    .section { background: white; border-radius: 12px; border: 1px solid #e8eaf0; margin-bottom: 1.25rem; overflow: hidden; }
    .section-header { padding: 0.75rem 1.5rem; border-bottom: 1px solid #e8eaf0; display: flex; align-items: center; gap: 8px; }
    .section-header h2 { font-size: 12px; font-weight: 600; color: #95a5a6; text-transform: uppercase; letter-spacing: 0.07em; }
    .section-icon { width: 16px; height: 16px; color: #95a5a6; flex: none; }
    .section-body { padding: 1.25rem 1.5rem; }
    .table-wrap { overflow-x: auto; }
    table { width: 100%; border-collapse: collapse; font-size: 13px; }
    th { text-align: left; padding: 8px 14px; font-size: 11px; font-weight: 600; color: #95a5a6; text-transform: uppercase; letter-spacing: 0.05em; border-bottom: 1px solid #e8eaf0; background: #fafbfc; white-space: nowrap; }
    th.grp { text-align: center; border-left: 1px solid #e8eaf0; }
    td { padding: 9px 14px; border-bottom: 1px solid #f0f2f5; color: #2c3e50; }
    th.num, td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
    tr:last-child td { border-bottom: none; }
    tbody tr:hover td { background: #fafbfc; }
    table.kv { width: auto; }
    table.kv th { background: none; border: none; padding: 4px 24px 4px 0; }
    table.kv td { border: none; padding: 4px 0; }
    .sub { color: #95a5a6; font-size: 12px; }
    .dash { color: #bdc3c7; }
    .lead { font-size: 13px; color: #52514e; margin-bottom: 0.75rem; }
    .note { font-size: 11px; color: #95a5a6; margin-top: 0.75rem; font-style: italic; }
    .chart { width: 100%; min-height: 380px; margin-bottom: 1rem; }
    .controls { display: flex; flex-wrap: wrap; align-items: flex-end; gap: 16px; margin-bottom: 1rem; }
    .control label { display: block; font-size: 11px; font-weight: 600; color: #95a5a6; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }
    .control select, .control input { font: inherit; font-size: 13px; color: #2c3e50; padding: 6px 10px; border: 1px solid #d5d8de; border-radius: 6px; background: white; }
    .control input { width: 110px; }
    .btn { font: inherit; font-size: 13px; padding: 7px 14px; border-radius: 6px; border: 1px solid #d5d8de; background: white; color: #2c3e50; cursor: pointer; }
    .btn:hover { background: #f5f6fa; }
    .btn-primary { background: #2980b9; border-color: #2980b9; color: white; }
    .btn-primary:hover { background: #2471a3; }
    .taxon-species { font-style: italic; }
    .swatch { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 8px; vertical-align: -1px; }
    tr.sample-start td { border-top: 1px solid #e8eaf0; }
    footer { background: white; border-top: 1px solid #e8eaf0; padding: 1rem 2.5rem; display: flex; flex-wrap: wrap; align-items: center; gap: 8px; }
    .footer-label { font-size: 11px; color: #95a5a6; margin-right: 4px; }
    .soft-tag { display: inline-block; padding: 2px 9px; background: #f8f9fa; border: 1px solid #e8eaf0; border-radius: 6px; font-size: 11px; color: #7f8c8d; }
    @media (max-width: 700px) { header, .container, footer { padding-left: 16px; padding-right: 16px; } }
</style>
</head>
<body>

<header>
    <div class="header-left">
        %%LOGO%%
        <div class="header-title">
            <h1>%%TITLE%%</h1>
            <p>Classification EMU · %%N_SAMPLES%% échantillons</p>
        </div>
    </div>
    <div class="header-right">
        <div class="date">%%DATE%%</div>
        <div class="version">%%VERSION%%</div>
    </div>
</header>

<div class="container">

    <div class="stats-grid">
        <div class="stat-card"><div class="value">%%N_SAMPLES%%</div><div class="label">Échantillons</div></div>
        <div class="stat-card"><div class="value">%%TOTAL_RAW%%</div><div class="label">Reads bruts</div></div>
        <div class="stat-card"><div class="value">%%TOTAL_FILT%% <small>%%PCT_FILT%%</small></div><div class="label">Reads après Chopper</div></div>
        <div class="stat-card"><div class="value">%%TOTAL_EMU%%</div><div class="label">Reads analysés par EMU</div></div>
    </div>

    %%PARAMS%%
    %%SUMMARY%%
    %%QC%%
    %%SPIKE%%

    <div class="section">
        <div class="section-header">
            <svg class="section-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 20V10m6 10V4m6 16v-7m6 7H2"/></svg>
            <h2>Abondances relatives (EMU)</h2>
        </div>
        <div class="section-body">
            <div class="controls">
                <div class="control">
                    <label for="rank">Rang</label>
                    <select id="rank">
                        <option value="species">Espèce</option>
                        <option value="genus">Genre</option>
                        <option value="family">Famille</option>
                    </select>
                </div>
                <div class="control">
                    <label for="minReads">Reads minimum</label>
                    <input id="minReads" type="number" min="0" step="1">
                </div>
                <div class="control">
                    <label for="minPct">Abondance minimale (%)</label>
                    <input id="minPct" type="number" min="0" max="100" step="0.1">
                </div>
                <button class="btn" id="resetBtn" type="button">Seuils par défaut</button>
                <button class="btn btn-primary" id="csvBtn" type="button">Télécharger le tableau (CSV)</button>
            </div>
            <div id="abundance_chart" class="chart"></div>
            <p class="note">Un taxon est affiché dans un échantillon s'il atteint <b>à la fois</b> le nombre de reads
                minimum <b>et</b> l'abondance minimale ; les autres sont regroupés dans « Sous les seuils ».
                Au-delà des 8 taxons colorés, les taxons au-dessus des seuils sont regroupés dans « Autres taxons »
                (détail au survol et dans le tableau). Abondances en % des reads classés par EMU, contrôle
                interne exclu.</p>
            <div class="table-wrap" id="abundance_table"></div>
        </div>
    </div>

</div>

<footer>
    <span class="footer-label">Outils :</span>
    %%TOOLS%%
    <span class="footer-label" style="margin-left:auto">Rapport généré le %%GENERATED%%</span>
</footer>

<script>
const QC = %%QC_DATA%%;
const AB = %%ABUNDANCE_DATA%%;

// Palette catégorielle validée (ordre fixe, jamais recyclée) + gris neutres
const PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"];
const OTHER_COLOR = "#898781";
const BELOW_COLOR = "#d5d4cd";
const N_COLORED = PALETTE.length;
const PLOT_CONFIG = {displaylogo: false, responsive: true,
                     modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"]};
const FONT = {family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif', size: 12, color: "#2c3e50"};

const fmtInt = v => Math.round(v).toLocaleString("fr-FR");
const fmtPct = v => v.toLocaleString("fr-FR", {minimumFractionDigits: 1, maximumFractionDigits: 1}) + " %";
const escHtml = s => String(s).replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));

// ── Graphique QC : reads avant / après Chopper (échelle log) ────────────────
Plotly.newPlot("qc_chart", [
    {type: "bar", name: "Avant Chopper", x: QC.samples, y: QC.raw, marker: {color: PALETTE[0], line: {color: "white", width: 1}},
     hovertemplate: "%{x}<br>Avant Chopper : %{y:,} reads<extra></extra>"},
    {type: "bar", name: "Après Chopper", x: QC.samples, y: QC.filtered, marker: {color: PALETTE[1], line: {color: "white", width: 1}},
     hovertemplate: "%{x}<br>Après Chopper : %{y:,} reads<extra></extra>"}
], {
    barmode: "group", bargap: 0.25, bargroupgap: 0.08, font: FONT, height: 380,
    margin: {t: 30, r: 10, b: 90, l: 70}, plot_bgcolor: "white", paper_bgcolor: "white",
    legend: {orientation: "h", x: 0, y: 1.12},
    yaxis: {type: "log", dtick: 1, exponentformat: "none", tickformat: ",d",
            title: {text: "Reads (échelle log)"}, gridcolor: "#ecf0f1", zeroline: false},
    xaxis: {tickangle: -45, automargin: true},
    separators: ", "
}, PLOT_CONFIG);

// ── Abondances ──────────────────────────────────────────────────────────────
const rankSel = document.getElementById("rank");
const minReadsIn = document.getElementById("minReads");
const minPctIn = document.getElementById("minPct");
minReadsIn.value = AB.minReads;
minPctIn.value = AB.minPct;

let currentRows = [];

function compute() {
    const rank = rankSel.value;
    const minReads = Math.max(0, parseFloat(minReadsIn.value) || 0);
    const minPct = Math.max(0, parseFloat(minPctIn.value) || 0);
    const taxa = AB.ranks[rank];
    const order = AB.order[rank];
    // Couleur liée au taxon (ordre fixe calculé une fois) : changer les
    // seuils ne repeint jamais un taxon déjà affiché.
    const colored = order.slice(0, N_COLORED);
    const colorOf = {};
    colored.forEach((t, i) => colorOf[t] = PALETTE[i]);

    const traces = colored.map(t => ({name: t, y: [], reads: []}));
    const others = {name: "Autres taxons", y: [], text: []};
    const below = {name: "Sous les seuils", y: [], text: []};
    const rows = [];

    AB.samples.forEach(s => {
        const total = AB.totals[s] || 0;
        const passing = [];
        let belowReads = 0, otherPct = 0;
        const otherNames = [];
        order.forEach(t => {
            const reads = (taxa[t] || {})[s] || 0;
            if (!reads || !total) return;
            const pct = 100 * reads / total;
            if (reads >= minReads && pct >= minPct) passing.push({t, reads, pct});
            else belowReads += reads;
        });
        traces.forEach(tr => {
            const hit = passing.find(p => p.t === tr.name);
            tr.y.push(hit ? hit.pct : 0);
            tr.reads.push(hit ? hit.reads : 0);
        });
        passing.filter(p => !(p.t in colorOf)).forEach(p => {
            otherPct += p.pct;
            otherNames.push(p.t + " : " + fmtPct(p.pct));
        });
        others.y.push(otherPct);
        others.text.push(otherNames.join("<br>"));
        const belowPct = total ? 100 * belowReads / total : 0;
        below.y.push(belowPct);
        below.text.push(fmtInt(belowReads) + " reads");

        passing.sort((a, b) => b.pct - a.pct);
        if (!total) {
            const label = (AB.classified[s] || 0) > 0 ? "Uniquement le contrôle interne" : "Aucun read classé";
            rows.push({sample: s, taxon: label, reads: null, pct: null, kind: "empty"});
        } else {
            passing.forEach(p => rows.push({sample: s, taxon: p.t, reads: p.reads, pct: p.pct, kind: "taxon"}));
            if (belowReads > 0) rows.push({sample: s, taxon: "Sous les seuils", reads: belowReads, pct: belowPct, kind: "below"});
        }
    });

    const data = traces.filter(tr => tr.y.some(v => v > 0)).map(tr => ({
        type: "bar", name: tr.name, x: AB.samples, y: tr.y, customdata: tr.reads,
        marker: {color: colorOf[tr.name], line: {color: "white", width: 1}},
        hovertemplate: "%{x}<br><b>" + escHtml(tr.name) + "</b> : %{y:.1f} % (%{customdata:,.0f} reads)<extra></extra>"
    }));
    if (others.y.some(v => v > 0)) data.push({
        type: "bar", name: others.name, x: AB.samples, y: others.y, text: others.text, textposition: "none",
        marker: {color: OTHER_COLOR, line: {color: "white", width: 1}},
        hovertemplate: "%{x}<br><b>Autres taxons</b> : %{y:.1f} %<br>%{text}<extra></extra>"
    });
    if (below.y.some(v => v > 0)) data.push({
        type: "bar", name: below.name, x: AB.samples, y: below.y, text: below.text, textposition: "none",
        marker: {color: BELOW_COLOR, line: {color: "white", width: 1}},
        hovertemplate: "%{x}<br><b>Sous les seuils</b> : %{y:.1f} % (%{text})<extra></extra>"
    });

    Plotly.react("abundance_chart", data, {
        barmode: "stack", bargap: 0.25, font: FONT, height: 460,
        margin: {t: 10, r: 10, b: 90, l: 60}, plot_bgcolor: "white", paper_bgcolor: "white",
        legend: {orientation: "v", x: 1.01, y: 1, font: {size: 11}},
        yaxis: {title: {text: "Abondance relative (%)"}, range: [0, 100], gridcolor: "#ecf0f1", zeroline: false},
        xaxis: {tickangle: -45, automargin: true},
        separators: ", "
    }, PLOT_CONFIG);

    renderTable(rows, rank, colorOf);
    currentRows = rows;
}

function renderTable(rows, rank, colorOf) {
    let html = "<table><thead><tr><th>Échantillon</th><th>Taxon</th><th class='num'>Reads estimés</th>" +
               "<th class='num'>Abondance</th></tr></thead><tbody>";
    let previous = null;
    rows.forEach(r => {
        const first = r.sample !== previous;
        previous = r.sample;
        let taxon = escHtml(r.taxon);
        if (r.kind === "taxon") {
            const color = colorOf[r.taxon] || OTHER_COLOR;
            taxon = "<span class='swatch' style='background:" + color + "'></span>" +
                    (rank === "species" ? "<span class='taxon-species'>" + taxon + "</span>" : taxon);
        } else {
            taxon = "<span class='sub'>" + taxon + "</span>";
        }
        html += "<tr" + (first ? " class='sample-start'" : "") + "><td>" + (first ? escHtml(r.sample) : "") + "</td>" +
                "<td>" + taxon + "</td>" +
                "<td class='num'>" + (r.reads === null ? "<span class='dash'>—</span>" : fmtInt(r.reads)) + "</td>" +
                "<td class='num'>" + (r.pct === null ? "<span class='dash'>—</span>" : fmtPct(r.pct)) + "</td></tr>";
    });
    document.getElementById("abundance_table").innerHTML = html + "</tbody></table>";
}

// Export CSV : séparateur « ; » et virgule décimale (Excel en français)
document.getElementById("csvBtn").addEventListener("click", () => {
    const num = (v, d) => v === null ? "" : v.toFixed(d).replace(".", ",");
    const quote = s => '"' + String(s).replace(/"/g, '""') + '"';
    const lines = ["Echantillon;Rang;Taxon;Reads estimes;Abondance (%)"];
    currentRows.forEach(r => lines.push([quote(r.sample), rankSel.value, quote(r.taxon), num(r.reads, 1), num(r.pct, 2)].join(";")));
    const blob = new Blob(["﻿" + lines.join("\r\n")], {type: "text/csv;charset=utf-8"});
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "abondances_" + rankSel.value + "_" + minReadsIn.value + "reads_" + minPctIn.value + "pct.csv";
    document.body.appendChild(link);
    link.click();
    link.remove();
});

document.getElementById("resetBtn").addEventListener("click", () => {
    minReadsIn.value = AB.minReads;
    minPctIn.value = AB.minPct;
    compute();
});
[rankSel, minReadsIn, minPctIn].forEach(el => el.addEventListener("input", compute));
compute();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
