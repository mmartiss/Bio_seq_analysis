from Bio import SeqIO

import argparse
import sys

# What needs to be done and how I visualize programs structure:
# 1. Input
# 2. Read 6 reading frames
# 3. Search for ORFs (start/stop codons) and their bounds
# 4. Evaluate and filter ORFs (What critieria?)
# 5. Print out results (Coordinates, sequence, reading frame, length, score)


# Testing bacteria genome
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_037044595.1/
# Mycoplasma mycoides JCVI-syn3.0
# Smallest genome


# Code recreated from:
# https://lanadominkovic.medium.com/bioinformatics-101-reading-fasta-files-using-biopython-501c390c6820
def fasta_in(file, debug):
    #for sequence in SeqIO.parse("input/JCVI-syn3.0.fasta", "fasta"):
    for sequence in SeqIO.parse(file, "fasta"):
        if debug:
            print(sequence.id)
            print(sequence.seq)
            print(repr(sequence.seq))
            print(len(sequence))
        seq = sequence.seq

    return str(seq)

def reverse_sequence(sequence):
    reverse_seq = []
    for n in list(sequence):
        if n == "A": reverse_seq.append("T")
        if n == "T": reverse_seq.append("A")
        if n == "C": reverse_seq.append("G")
        if n == "G": reverse_seq.append("C")

    reverse_seq.reverse()
    return "".join(reverse_seq)

def reading_frame(start, sequence):
    if start+3 <= len(sequence) and start >= 0:
        frame = sequence[start:start+3]
    else:
        return None
    return frame

def reading_frame_shift(shift, sequence):
    # iteruoti per butent ta shift
    shifted_frame = []
    for i in range(shift, len(sequence), 3):
        tmp_frame=[]
        tmp_frame.append(reading_frame(i, sequence))
        tmp_frame.append(shift)
        tmp_frame.append(i)
        shifted_frame.append(tmp_frame)

    return shifted_frame


def get_all_possible_reading_frames(sequence):
    # reikia poslinkiu 0,1,2
    # reikia indekso keitimo kas 3
    tmp_frames_list=[]
    frames=[]
    for shift in range(3):
        tmp_frames_list = reading_frame_shift(shift, sequence)
        frames = frames + tmp_frames_list

    return frames

def remove_none_from_list(frame_list):
    new_list = list(filter(lambda x: x[0] is not None, frame_list))
    return new_list

def find_possible_orfs(frames, strand, sequence, total_genome_len):
    which = 0
    possible_orfs = []
    for frame in frames:
        if frame[0] == 'ATG' or frame[0] == 'GTG' or frame[0] == 'TTG':
            for stop in frames:
                if frame[1] == stop[1] and frame[2] < stop[2] and (stop[0] == 'TAA' or stop[0] == 'TAG' or stop[0] == 'TGA'):
                    seq = sequence[frame[2]:stop[2]+3]
                    length = len(seq)//3 # nes antraip float

                    if strand == "+":
                        g_start = frame[2]
                        g_end = stop[2] + 3
                    else:  # "-"
                        g_start = total_genome_len - (stop[2] + 3)
                        g_end = total_genome_len - frame[2]

                    orf = {
                        "start": frame[0],
                        "stop": stop[0],
                        "shift": frame[1],
                        "strand": strand,
                        "genomic_start": g_start,
                        "genomic_end": g_end,
                        "sequence": seq,
                        "length": length,
                        "index" : which
                    }
                    possible_orfs.append(orf)
                    which += 1

    return possible_orfs

def filter_overlapping_orfs(orfs_list, max_allowed_overlap=15):
    sorted_by_length = sorted(orfs_list, key=lambda x: x["length"], reverse=True)
    
    clean_orfs = []

    for candidate in sorted_by_length:
        keep_candidate = True

        for accepted in clean_orfs:
            if candidate["genomic_start"] < accepted["genomic_end"] and accepted["genomic_start"] < candidate["genomic_end"]:
                overlap_len = min(candidate["genomic_end"], accepted["genomic_end"]) - max(candidate["genomic_start"], accepted["genomic_start"])
                
                if overlap_len > max_allowed_overlap:
                    keep_candidate = False
                    break
        
        if keep_candidate:
            clean_orfs.append(candidate)
            
    return sorted(clean_orfs, key=lambda x: x["genomic_start"])

# https://pmc.ncbi.nlm.nih.gov/articles/PMC3454372/
def filter_wrong_length_orfs(orfs_list):
    correct_orfs = []
    for orf in orfs_list:
        if 100 <= orf["length"] <= 700:
            correct_orfs.append(orf)
            
    return correct_orfs

def filtering_orfs(orfs_list):
    correct_orfs = filter_wrong_length_orfs(orfs_list)
    correct_orfs = filter_overlapping_orfs(correct_orfs, max_allowed_overlap=15)

    return correct_orfs

def main():
    parser = argparse.ArgumentParser(
        description="Find possible ORFs"
    )
    parser.add_argument(
        "-i", "--input", type=str, help="FIle name/path"
    )
    parser.add_argument(
        "-d", "--default-file", action="store_true",
        help="Run in test mode using input/test.fasta",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable debug mode"
    )

    args = parser.parse_args()
    file = args.input

    if not file:
        if args.default_file:
            file = "input/test.fasta"
        else:
            print("No input file specified via -i / --input.")
            sys.exit(1)

    print(f"Selected file: {file}")
    print(f"Debug mode: {'ENABLED' if args.debug else 'DISABLED'}")

    #Original sequence
    input_seq = fasta_in(file, args.debug)
    N = len(input_seq)
    reverse_seq = reverse_sequence(input_seq)
    forward = get_all_possible_reading_frames(input_seq)
    reverse = get_all_possible_reading_frames(reverse_seq)
    forward_clean_frames = remove_none_from_list(forward)
    reverse_clean_frames = remove_none_from_list(reverse)

    if args.debug:
        print("\nAll frames\n")
        print("Forward cleaned frames")
        print(forward_clean_frames)
        print("\nReverse cleaned frames")
        print(reverse_clean_frames)

    f_orfs = find_possible_orfs(forward_clean_frames, "+", input_seq, N)
    r_orfs = find_possible_orfs(reverse_clean_frames, "-", reverse_seq, N)

    if args.debug:
        print("\nAll possible frames \n")
        print("Forward")
        print(f_orfs)
        print("\nReverse")
        print(r_orfs)

    all_orfs = f_orfs + r_orfs

    good_orfs = filtering_orfs(all_orfs)

    print("---------------------------------------------------------------")
    print("                             GOOD:D                             ")
    print("---------------------------------------------------------------")
    print(good_orfs)

if __name__ == '__main__':
    main()