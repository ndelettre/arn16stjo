#!/usr/bin/env python3
"""
subsample_fastq.py
------------------
Tire au hasard au plus N reads d'un FASTQ (gzip ou non), de façon
reproductible (graine fixe) : relancer le pipeline donne les mêmes reads.

- Si le FASTQ contient N reads ou moins, tous les reads sont gardés.
- Sinon, N reads sont tirés sans remise ; leur ordre d'origine est conservé.
- N = 0 désactive le sous-échantillonnage (tous les reads sont gardés).
- Sans --output, les reads sont seulement comptés (aucune copie écrite).

Deux lectures du fichier : la 1re compte les reads, la 2de écrit ceux qui ont
été tirés. La mémoire utilisée reste faible, même pour des millions de reads.

Écrit aussi un TSV : sample | reads_nohost | reads_emu | subsampled

Usage :
    subsample_fastq.py --input s_nohost.fastq.gz --max-reads 20000 --seed 42 \
        --sample s --output s_emu_input.fastq.gz --stats s_subsample.tsv
"""

import argparse
import gzip
import random


def open_fastq(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def count_reads(path):
    with open_fastq(path) as handle:
        return sum(1 for _ in handle) // 4


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="FASTQ d'entrée (.fastq ou .fastq.gz)")
    parser.add_argument("--max-reads", required=True, type=int, help="nombre maximal de reads (0 = tous)")
    parser.add_argument("--seed", type=int, default=42, help="graine du tirage aléatoire [42]")
    parser.add_argument("--sample", required=True, help="nom de l'échantillon")
    parser.add_argument("--output", help="FASTQ de sortie (.fastq, ou .fastq.gz pour compresser) ; "
                                         "sans --output, les reads sont seulement comptés")
    parser.add_argument("--stats", required=True, help="TSV récapitulatif")
    args = parser.parse_args()

    n_reads = count_reads(args.input)
    subsample = 0 < args.max_reads < n_reads
    if subsample:
        keep = set(random.Random(args.seed).sample(range(n_reads), args.max_reads))
    n_kept = args.max_reads if subsample else n_reads

    if args.output:
        if args.output.endswith(".gz"):
            fout_handle = gzip.open(args.output, "wt", compresslevel=1)
        else:
            fout_handle = open(args.output, "w")
        with open_fastq(args.input) as fin, fout_handle as fout:
            for index in range(n_reads):
                record = [next(fin) for _ in range(4)]
                if not subsample or index in keep:
                    fout.writelines(record)

    with open(args.stats, "w", encoding="utf-8") as out:
        out.write("sample\treads_nohost\treads_emu\tsubsampled\n")
        out.write(f"{args.sample}\t{n_reads}\t{n_kept}\t{'oui' if subsample else 'non'}\n")


if __name__ == "__main__":
    main()
