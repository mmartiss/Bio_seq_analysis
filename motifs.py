from Bio import SeqIO
from collections import Counter

import argparse
import sys

def fasta_in(file, debug):
    seq = []
    for sequences in SeqIO.parse(file, "fasta"):
        if debug:
            print(sequences.id)
            print(sequences.seq)
            print(repr(sequences.seq))
            print(len(sequences))
        seq.append(str(sequences.seq))

    return seq

def find_k_mers(seqs, k):
    k_mers_list = []

    for seq in seqs:
        k_mers = []
        for i in range(len(seq) - k + 1):
            possible_k_mer = seq[i:i + k]
            k_mers.append(str(possible_k_mer))

        k_mers_list.append(k_mers)

    return k_mers_list

def k_mer_count(k_mers, debug):
    total_counts = Counter()
    sequence_counts = []

    for index, sequence_k_mers in enumerate(k_mers, start=1):
        counts = Counter(sequence_k_mers)
        sequence_counts.append(counts)
        total_counts.update(counts)

        if debug:
            print("Sequence: ", index)
            for k_mer, count in counts.items():
                print(k_mer, count)

            print("Total: ", "\n")
            for k_mer, count in total_counts.items():
                print(k_mer, count)

    return total_counts, sequence_counts

def main():
    parser = argparse.ArgumentParser(
        description="Find motifs"
    )
    parser.add_argument(
        "-i", "--input", type=str, help="File name/path"
    )
    parser.add_argument(
        "-d", "--default-file", action="store_true",
        help="Run in test mode using input/test.fasta",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable debug mode"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output.fasta", help="Output file name or path"
    )

    args = parser.parse_args()
    file = args.input

    if not file:
        if args.default_file:
            file = "input/test_motifs.fasta"
        else:
            print("No input file specified via -i / --input.")
            sys.exit(1)

    print(f"Selected file: {file}")
    print(f"Debug mode: {'ENABLED' if args.debug else 'DISABLED'}")

    input_seqs = fasta_in(file, args.debug)
    if args.debug:
        print(input_seqs)

    k = 3
    k_mers = find_k_mers(input_seqs, k)
    if args.debug:
        print(k_mers)

    total_counts, k_mer_counts = k_mer_count(k_mers, args.debug)
    if args.debug:
        print(k_mer_counts)

if __name__ == '__main__':
    main()