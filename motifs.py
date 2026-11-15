from Bio import SeqIO
from collections import Counter

import argparse
import sys
import random

# INPUT: file - FASTA file path, debug - true/false
# OUTPUT: seq - list of nucleotide sequence strings
# fasta_in function reads all sequences from a FASTA file using SeqIO
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

# INPUT: seqs - list of sequences, k - k-mer length
# OUTPUT: k_mers_list - list of k-mer lists, one list per sequence
# find_k_mers function extracts all overlapping k-mers while keeping their sequence of origin
def find_k_mers(seqs, k):
    k_mers_list = []

    for seq in seqs:
        k_mers = []
        for i in range(len(seq) - k + 1):
            possible_k_mer = seq[i:i + k]
            k_mers.append(str(possible_k_mer))

        k_mers_list.append(k_mers)

    return k_mers_list

# INPUT: k_mers - list of k-mer lists, debug - true/false
# OUTPUT: total_counts - dictionary of total counts, 
# sequence_counts - list of per-sequence count dictionaries
# k_mer_count function counts k-mers within each sequence and across all sequences
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


# INPUT: input_seqs - iterable of nucleotide strings, debug - true/false
# OUTPUT: n_counts - dictionary of A, T, C and G counts
# nucleotide_counts function counts uppercase A, T, C and G across the supplied strings
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

# INPUT: a, b, c, d - nucleotide counts with a positive total
# OUTPUT: relative frequency of a (float)
# frequencies function divides the first count by the sum of all four counts
def frequencies(a, b, c, d):
    return a/(a+b+c+d)

# INPUT: n_counts - dictionary of A, T, C and G counts, debug - true/false
# OUTPUT: n_freqs - dictionary of nucleotide frequencies
# nucleotide_frequencies function converts nucleotide counts into frequencies using their total
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

# INPUT: n_freqs - background nucleotide frequencies, k_mer - nucleotide string, debug - true/false
# OUTPUT: prob - background probability of the k-mer (float)
# k_mer_prob function multiplies background probabilities assuming independent nucleotides
def k_mer_prob(n_freqs, k_mer, debug):
    prob = 1
    for n in k_mer:
        prob = prob * n_freqs[n]

    if debug:
        print(k_mer, " : ", prob) 

    return prob

# INPUT: n_freqs - background frequencies, k_mers - list of k-mer lists, debug - true/false
# OUTPUT: p - dictionary mapping each observed distinct k-mer to its background probability
# probabilities function calculates background probabilities for the observed k-mers
def probabilities(n_freqs, k_mers, debug):
    p = {}

    for m in k_mers:
        for mer in m:
            #p.append(k_mer_prob(n_freqs, mer, debug))
            p[mer] = k_mer_prob(n_freqs, mer, debug)

    if debug:
        print("Probabilities: \n", p)

    return p

# INPUT: k_mers - list of k-mer lists, debug - true/false
# OUTPUT: total_windows - total number of k-mer windows (integer)
# count_k_mer_windows function adds the number of windows from every sequence
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

# INPUT: total_counts - observed k-mer counts, prob - background probabilities, 
# total_win - total windows, debug - true/false
# OUTPUT: rez - dictionary of observed counts, expected counts and enrichment values
# k_mer_enrich function calculates expected counts and observed/expected ratios 
# uses None when expected is zero
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

# INPUT: item - k-mer and result dictionary pair, debug - true/false
# OUTPUT: enrichment value, or -1 when enrichment is None
# get_enrich function extracts an enrichment value for sorting results
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
# INPUT: results - enrichment results, sequence_counts - per-sequence counts, 
# debug - true/false, output - output file path
# OUTPUT: STDOUT and file
# print_results function sorts k-mers by descending enrichment and saves the table, 
# overwriting existing file contents
def print_results(results, sequence_counts, debug, output):

    # I'm sorry for this nested function, just wanted to remember that
    # INPUT: item - k-mer and result dictionary pair. debug comes from the enclosing function
    # OUTPUT: sorting value returned by get_enrich
    # sorting_key function provides the enrichment sorting key for the results table
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


# INPUT: seqs - list of sequences, k - motif length, debug - true/false
# OUTPUT: pos - list of zero-based motif starting positions
# rand_pos function chooses one valid random motif starting position in each sequence
def rand_pos(seqs, k, debug):
    # Raise ValueError, wanted to do it earlier, but jsut got to it
    # Moved to main
    pos = []

    if debug:
        print("Random positions:")
    for seq in seqs:
        p = (random.randint(0, len(seq) - k))
        if debug:
            print(seq, " : ", p, ":", seq[p:p + k])

        pos.append(p)

    return pos

# INPUT: seqs - list of sequences, debug - true/false
# OUTPUT: selected - randomly selected zero-based sequence index
# rand_seq function chooses the sequence to leave out when building the motif profile
def rand_seq(seqs, debug):
    selected = random.randrange(len(seqs))

    if debug:
        print(f"Selected random seq: {selected}")

    return selected

# INPUT: seqs - sequences, selected - index to exclude, pos - motif positions, 
# k - motif length, debug - true/false
# OUTPUT: motifs - list of current motif strings excluding the selected sequence
# get_motifs function extracts the current motifs from all sequences except the selected one
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

# INPUT: motifs - equal-length motif strings, k - motif length, debug - true/false
# OUTPUT: counts - list of nucleotide count dictionaries, one per motif position
# motif_counts function counts A, T, C and G separately in each motif column
def motif_counts(motifs, k, debug):
    counts = []

    for position in range(k):
        column = []

        for motif in motifs:
            column.append(motif[position])

        # Karma points for reusing code:D
        column_counts = nucleotide_counts(column, False)
        counts.append(column_counts)

    if debug:
        for position, column_counts in enumerate(counts):
            print(f"Position {position}: {column_counts}")

    return counts

# INPUT: counts_matrix - per-position nucleotide counts, debug - true/false
# OUTPUT: prof - list of per-position nucleotide probability dictionaries
# motifs_prof function adds a pseudocount of one to each nucleotide count and normalizes each column
def motifs_prof(counts_matrix, debug):
    prof = []

    for column_counts in counts_matrix:
        adjusted_c = {} # I think it sould be a dict
        # easier to access ig
        for nucl, count in column_counts.items():
            adjusted_c[nucl] = count + 1

        # Haha more reusable code:D
        column_freqs = nucleotide_frequencies(adjusted_c, debug)
        prof.append(column_freqs)

    if debug:
        for pos, freqs in enumerate(prof):
            print(f"Position {pos}: {freqs}")

    return prof

# INPUT: seq - held-out sequence, k - motif length, profile - positional probabilities, 
# background - nucleotide frequencies, debug - true/false
# OUTPUT: weights - list of sampling weights indexed by candidate start position
# motif_weights function scores every candidate window using the product of 
# motif/background probability ratios
def motif_weights(seq, k, profile, background, debug):
    weights = []

    for start in range(len(seq) - k + 1):
        motif = seq[start:start + k]
        weight = 1 # nto zero

        for position, nucleotide in enumerate(motif):
            weight *= (profile[position][nucleotide] / background[nucleotide])

        weights.append(weight)

        if debug:
            print(f"Position {start}: motif={motif}, weight={weight}")

    return weights

# INPUT: weights - sampling weights for all candidate positions, debug - true/false
# OUTPUT: selected - sampled zero-based motif starting position
# sample_pos function randomly samples a position in proportion to its weight
def sample_pos(weights, debug):
    pos = range(len(weights))

    # I didnt knwo how to do this better or more understandable
    # Bc we need to choose the better weight choise, but
    # still need a way to be able to choose a different weight
    # which may not bet better
    # Ig this is how to leave a posibility to try other positionss
    selected = random.choices(pos, weights=weights, k=1)[0]

    if debug:
        print("New motif pos: ", selected)

    return selected

# AI disclaimer, wrote this code mysel, but asked AI to make it into a function.
# Didn't want it in the main function as masive code snippet.
# So AI rewrote variable names and made this into a function.
# INPUT: seqs - sequences, freqs - background frequencies, k - motif length, 
# iterations - update count, debug - true/false
# OUTPUT: best_positions - best position list, best_score - lowest mismatch score found
# pos_loop function runs one Gibbs search and retains the initial 
# or sampled solution with the fewest consensus mismatches.
def pos_loop(seqs, freqs, k, iterations, debug):
    positions = rand_pos(seqs, k, debug)

    best_positions = positions.copy()
    best_score = motif_score(seqs, positions, k, debug)

    for iteration in range(iterations):
        if debug:
            print(f"\nIteration {iteration + 1}")

        rand_sequence = rand_seq(seqs, debug)

        motifs = get_motifs(
            seqs, rand_sequence, positions, k, debug
        )

        counts_matrix = motif_counts(motifs, k, debug)
        profile = motifs_prof(counts_matrix, debug)

        weights = motif_weights(
            seqs[rand_sequence],
            k,
            profile,
            freqs,
            debug
        )

        positions[rand_sequence] = sample_pos(weights, debug)

        score = motif_score(seqs, positions, k, False)

        if score < best_score:
            best_score = score
            best_positions = positions.copy()

            if debug:
                print(f"New best score: {best_score}")

    return best_positions, best_score

# INPUT: seqs - sequences, positions - motif starts, k - motif length, debug - true/false
# OUTPUT: score - total number of mismatches to the column consensus (integer)
# motif_score function evaluates motifs from all sequences 
# lower scores mean fewer consensus mismatches.
def motif_score(seqs, positions, k, debug):
    motifs = []

    for seq, position in zip(seqs, positions):
        motifs.append(seq[position:position + k])

    counts_matrix = motif_counts(motifs, k, debug)

    score = 0

    for column_counts in counts_matrix:
        most_common = max(column_counts.values())
        score += len(motifs) - most_common

        if debug:
            print("Most common: ", most_common)
            print("Score: ", score)

    return score

# AI disclaimer
# This function was written with the help of ai.
# It suggested this programmin decision for me
# To use inf for restarts for gibbs
# AI suggestion: use float("inf"), so first run score would be lower and it would become best starting result
# INPUT: seqs - sequences, freqs - background frequencies, k - motif length, 
# iterations - updates per run, restarts - number of runs, debug - true/false
# OUTPUT: best_positions - best position list across runs, best_score - lowest mismatch score
# gibbs_search function repeats Gibbs search from new random positions and keeps the best result
def gibbs_search(seqs, freqs, k, iterations, restarts, debug):
    best_positions = None
    best_score = float("inf")

    for restart in range(restarts):
        positions, score = pos_loop(
            seqs, freqs, k, iterations, debug
        )

        if debug:
            print(f"Restart {restart + 1}: best score={score}")

        if score < best_score:
            best_score = score
            best_positions = positions.copy()

    return best_positions, best_score

# AI disclaimer
# Asked AI to make output pretier
# INPUT: seqs - sequences, pos - best motif starts, k - motif length, 
# score - best mismatch score, output - file path
# OUTPUT: STDOUT, prints the Gibbs table and appends it to the output file
# gibbs_results function formats motif sequences and 
# zero-based positions, then appends them after the first results table
def gibbs_results(seqs, pos, k, score, output):
    lines = [
        "",
        "--------------------------------------------------------",
        "Gibbs sampling",
        f"Best motif score: {score}",
        "Positions are zero-based",
        "--------------------------------------------------------",
        f"{'Sequence':<12}{'Position':>10}  {'Motif'}"
    ]

    for i, seq in enumerate(seqs):
        position = pos[i]
        motif = seq[position:position + k]

        lines.append(
            f"{i + 1:<12}{position:>10}  {motif}"
        )

    table = "\n".join(lines)

    print(table)

    with open(output, "a", encoding="utf-8") as output_file:
        output_file.write(table + "\n")
# -------------------------------------------------------------------------------------------------
# INPUT: STDIN and file
# OUTPUT: STDOUT
# main function parses options, seeds randomness, 
# validates sequences and runs background enrichment and Gibbs motif searches
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
        "-o", "--output", type=str, default="motifs_output.txt", help="Output file name or path"
    )
    parser.add_argument(
        "-k", "--kmer",  type=int, default=3, help="k count"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducible results"
    )
    parser.add_argument(
        "--iter", type=int, default=1000,
        help="Number of gibbs sampling iterations"
    )
    parser.add_argument(
        "--restarts", type=int, default=10,
        help="Number of independent gibbs sampling runs"
    )

    args = parser.parse_args()
    file = args.input
    random.seed(args.seed)

    if args.iter < 1 or args.restarts < 1:
        parser.error("--iter and --restarts must be positive")

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
    print(f"Gibbs sampling iterations: {args.iter}, default = 1000")
    print("--------------------------------------------------------")
    print("Extra INFO:")
    print("If you wnat to change or havent changed the output file")
    print("you can do it wiht -o or --output flags")
    print("motifs_output.fasta is the default output file.")
    print("--------------------------------------------------------")

    input_seqs = fasta_in(file, args.debug)
    if args.debug:
        print(input_seqs)

    # -------------------------------------------------------------------------------------------------
    # Extra checks
    if len(input_seqs) < 2:
        parser.error("At least two sequences are required")

    if k < 1 or any(len(seq) < k for seq in input_seqs):
        parser.error("k must be positive and fit every sequence")

    if any(set(seq) - set("ACGT") for seq in input_seqs):
        parser.error("Sequences must contain only A, C, G and T")
    # -------------------------------------------------------------------------------------------------

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

    positions, score = gibbs_search(
        input_seqs, freqs, k, args.iter, args.restarts, args.debug
    )

    gibbs_results(input_seqs, positions, k, score, args.output)
    print("--------------------------------------------------------")


    # positions = rand_pos(input_seqs, k, args.debug)
    # rand_sequence = rand_seq(input_seqs, args.debug) # seq to exclude
    # motifs = get_motifs(input_seqs, rand_sequence, positions, k, args.debug)

    # counts_matrix = motif_counts(motifs, k, args.debug)

    # if args.debug:
    #     print("Motif lengths:", [len(motif) for motif in motifs])
    #     print("Matrix columns:", len(counts_matrix))

    # # profile[0]["A"] - A probability in first position of the motif
    # # IF you can call it like that
    # profile = motifs_prof(counts_matrix, args.debug)

    # weights = motif_weights(input_seqs[rand_sequence], k, profile, freqs, args.debug)

    # new_pos = pos_loop(input_seqs, freqs, k, args.iter, args.debug)

    # Commented out for dublication and test purposes.
    # positions, score = pos_loop(
    #     input_seqs, freqs, k, args.iter, args.debug
    # )

    # print(f"\nBest motif score: {score}")
    # gibbs_results(input_seqs, positions, k, args.debug)

    # score = motif_score(input_seqs, positions, k, args.debug)
    # if args.debug:
    #     print(score)

if __name__ == '__main__':
    main()