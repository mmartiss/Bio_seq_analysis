from Bio import SeqIO
from collections import Counter

import argparse
import sys
import random

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

# -------------------------------------------------------------------------------------------------

def find_k_mers(seqs, k):
    k_mers_list = []

    for seq in seqs:
        k_mers = []
        for i in range(len(seq) - k + 1):
            possible_k_mer = seq[i:i + k]
            k_mers.append(str(possible_k_mer))

        k_mers_list.append(k_mers)

    return k_mers_list

def k_mer_count(k_mers, debug=False):
    total_counts = Counter()
    sequence_counts = []

    for index, sequence_k_mers in enumerate(k_mers, start=1):
        counts = Counter(sequence_k_mers)
        sequence_counts.append(dict(counts))
        total_counts.update(counts)

        if debug:
            print("Sequence:", index)
            for k_mer, count in counts.items():
                print(k_mer, count)

    if debug:
        print("\nTotal:")
        for k_mer, count in total_counts.items():
            print(k_mer, count)

    return dict(total_counts), sequence_counts

# -------------------------------------------------------------------------------------------------


def nucleotide_counts(input_seqs, debug):
    a = t = c = g = 0

    for seq in input_seqs:
        for n in seq:
            if n == "A": a+=1
            if n == "T": t+=1
            if n == "C": c+=1
            if n == "G": g+=1

    n_counts = {"A": a, "T": t, "C": c, "G": g}
    if debug:
        print(n_counts)
    return n_counts

def frequencies(a, b, c, d):
    return a/(a+b+c+d)

def nucleotide_frequencies(n_counts, debug):
    n_freqs = {
        "A": frequencies(n_counts["A"], n_counts["T"], n_counts["C"], n_counts["G"]),
        "T": frequencies(n_counts["T"], n_counts["A"], n_counts["C"], n_counts["G"]),
        "C": frequencies(n_counts["C"], n_counts["A"], n_counts["T"], n_counts["G"]),
        "G": frequencies(n_counts["G"], n_counts["A"], n_counts["T"], n_counts["C"])
        }

    if debug:
        print(n_freqs)
    return n_freqs

def k_mer_prob(n_freqs, k_mer, debug):
    prob = 1
    for n in k_mer:
        prob = prob * n_freqs[n]

    if debug:
        print(k_mer, " : ", prob) 

    return prob

def probabilities(n_freqs, k_mers, debug):
    p = {}

    for m in k_mers:
        for mer in m:
            #p.append(k_mer_prob(n_freqs, mer, debug))
            p[mer] = k_mer_prob(n_freqs, mer, debug)

    if debug:
        print("Probabilities: \n", p)

    return p

def count_k_mer_windows(k_mers, debug):
    total_windows = 0

    for index, sequence_k_mers in enumerate(k_mers, start=1):
        sequence_windows = len(sequence_k_mers)
        total_windows += sequence_windows

        if debug:
            print(f"Sequence {index}: {sequence_windows} windows")

    if debug:
        print("Total windows: ", total_windows)

    return total_windows

def k_mer_enrich(total_counts, prob, total_win, debug):
    rez = {}

    for k_mer, observed in total_counts.items():
        expect = total_win * prob[k_mer]
        enrich = observed / expect if expect > 0 else None

        rez[k_mer] = {
            "observed": observed,
            "expected": expect,
            "enrichment": enrich
        }

        if debug:
            print(f"{k_mer}: o={observed}, ex={expect:.3f}, en={enrich}")
            #print(k_mer, " : ", observed, " : ", expect.3f, " : ", enrich, "\n")

    return rez

def get_enrich(item, debug):
    value = item[1]
    enrich = value["enrichment"]

    if enrich is None:
        if debug:
            print(f"{value}: changed None to -1\n")
        return -1

    return enrich

# Disclaimer, this function was written with the help of AI
# I didn't know how to make the output look GOOD.
# So AI helped with that.
def print_results(results, sequence_counts, debug, output):

    # I'm sorry for this nested function, just wanted to remember that
    def sorting_key(item):
        return get_enrich(item, debug)

    sorted_rez = sorted(
        results.items(),
        key=sorting_key,
        # I think it could look like this:
        # key=lambda item: get_enrich(item, debug),
        reverse=True
    )

    #print("\nK-mer\tObserved\tExpected\tEnrichment\tSequences")

    lines = [
        f"{'K-mer':<10}"
        f"{'Observed':>10}"
        f"{'Expected':>12}"
        f"{'Enrichment':>12}"
        f"{'Sequences':>12}"
    ]

    for k_mer, values in sorted_rez:
        sequence_number = sum(
            1 for counts in sequence_counts if k_mer in counts
        )

        enrichment = values["enrichment"]
        enrichment_text = (
            f"{enrichment:.3f}" if enrichment is not None else "NA"
        )

        sequence_text = f"{sequence_number}/{len(sequence_counts)}"

        line = (
            f"{k_mer:<10}"
            f"{values['observed']:>10}"
            f"{values['expected']:>12.3f}"
            f"{enrichment_text:>12}"
            f"{sequence_text:>12}"
        )

        lines.append(line)

    table = "\n".join(lines)

    print("--------------------------------------------------------")
    print("Motifs looked up with background")
    print("--------------------------------------------------------")
    print(table)
    print("--------------------------------------------------------")

    with open(output, "w", encoding="utf-8") as output_file:
        output_file.write(table + "\n")

# -------------------------------------------------------------------------------------------------

# Gibbs sampling

# https://www.geeksforgeeks.org/machine-learning/what-is-gibbs-sampling/
# https://www.youtube.com/watch?v=NRjhfyXWHuQ

# Raise ValueError, wanted to do it earlier, but jsut got to it
def rand_pos(seqs, k, debug):
    pos = []

    if debug:
        print("Random positions:")
    for seq in seqs:
        if not 1<= k <= len(seq):
            raise ValueError("k must be postitive and fit every seq")

        p = (random.randint(0, len(seq) - k))
        if debug:
            print(seq, " : ", p, ":", seq[p:p + k])

        pos.append(p)

    return pos

def rand_seq(seqs, debug):
    selected = random.randrange(len(seqs))

    if debug:
        print(f"Selected random seq: {selected}")

    return selected

def get_motifs(seqs, selected, pos, k, debug):
    motifs = []

    for i, seq in enumerate(seqs):
        if i == selected:
            if debug:
                print(f"Selected sequence to exlude: {i}, {seq}")
            continue

        position = pos[i]
        motif = seq[position:position + k]
        motifs.append(motif)

        if debug:
            print(f"Index: {i}, sequence: {seq}, pos: {position}, motif: {motif}")

    return motifs

# Karma points for reusing code:D
def motif_counts(motifs, k, debug):
    counts = []

    for position in range(k):
        column = []

        for motif in motifs:
            column.append(motif[position])

        column_counts = nucleotide_counts(column, False)
        counts.append(column_counts)

    if debug:
        for position, column_counts in enumerate(counts):
            print(f"Position {position}: {column_counts}")

    return counts

# -------------------------------------------------------------------------------------------------
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
        "-o", "--output", type=str, default="motifs_output.fasta", help="Output file name or path"
    )
    parser.add_argument(
        "-k", "--kmer",  type=int, default=3, help="k count"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducible results"
    )

    args = parser.parse_args()
    file = args.input
    random.seed(args.seed)

    if not file:
        if args.default_file:
            file = "input/test_motifs.fasta"
        else:
            print("No input file specified via -i / --input.")
            sys.exit(1)

    k = args.kmer

    print("--------------------------------------------------------")
    print(f"Selected file: {file}")
    print(f"Output file: {args.output}")
    print(f"Debug mode: {'ENABLED' if args.debug else 'DISABLED'}")
    print(f"K count: {k}")
    print(f"Seed: {args.seed}")
    print("--------------------------------------------------------")
    print("Extra INFO:")
    print("If you wnat to change or havent changed the output file")
    print("you can do it wiht -o or --output flags")
    print("motifs_output.fasta is the default output file.")
    print("--------------------------------------------------------")

    input_seqs = fasta_in(file, args.debug)
    if args.debug:
        print(input_seqs)

    k_mers = find_k_mers(input_seqs, k)
    if args.debug:
        print(k_mers)

    total_counts, k_mer_counts = k_mer_count(k_mers, args.debug)

    if args.debug:
        print("Seperate sequences counts: \n", k_mer_counts, "\n")
        print("All counts: \n", total_counts, "\n")

    counts = nucleotide_counts(input_seqs, args.debug)
    freqs = nucleotide_frequencies(counts, args.debug)

    k_mer_prob(freqs, k_mers[0][0], args.debug)

    prob = probabilities(freqs, k_mers, args.debug)

    total_windows = count_k_mer_windows(k_mers, args.debug)

    results = k_mer_enrich(total_counts, prob, total_windows, args.debug)

    print_results(results, k_mer_counts, args.debug, args.output)

    # Background end

    # Gibbs start:

    positions = rand_pos(input_seqs, k, args.debug)
    rand_sequence = rand_seq(input_seqs, args.debug) # seq to exclude
    motifs = get_motifs(input_seqs, rand_sequence, positions, k, args.debug)
    counts_matrix = motif_counts(motifs, k, args.debug)

if __name__ == '__main__':
    main()