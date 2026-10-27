from Bio import SeqIO

import argparse
import sys

# USED CRITERIA FOR POSSIBLE ORFs
# Length : from 100 to 700 of codons
# Overlapping: max 15
# Start codons: ATG, GTG, TTG
# Stop codons: TAA, TGA, TAG

# Testing bacteria genome
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_037044595.1/
# Mycoplasma mycoides JCVI-syn3.0
# Smallest genome
# Breaks my pc :D


# Code recreated from (input):
# https://lanadominkovic.medium.com/bioinformatics-101-reading-fasta-files-using-biopython-501c390c6820

# INPUT: file - file path, debug - true/false
# OUTPUT: seq - fasta sequence read from the file (string)
# fasta_in function reads fasta file using SeqIO parser. RETURNS ONLY THE LAST SEQUENCE!
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

# INPUT: sequence - nucleotide sequence
# OUTPUT: reverse_seq - created complementary reverse sequence
# reverse_sequence function creates complementary reverse sequence from original sequence
def reverse_sequence(sequence):
    reverse_seq = []
    for n in list(sequence):
        if n == "A": reverse_seq.append("T")
        if n == "T": reverse_seq.append("A")
        if n == "C": reverse_seq.append("G")
        if n == "G": reverse_seq.append("C")

    reverse_seq.reverse()
    return "".join(reverse_seq)


# INPUT: start - position, sequence
# OUTPUT: frame - a 3 nucleotide frame
# reading_frame gets a reading frame of 3 if possible, otherwise returns None
def reading_frame(start, sequence):
    if start+3 <= len(sequence) and start >= 0:
        frame = sequence[start:start+3]
    else:
        return None
    return frame

# INPUT: shift - 0 or 1 or 2 (reading frames shift), sequence
# OUPTUT: shifter_frame - all reading frames for that shift, its a list with info on seq, shift and position
# reading_frame_shift gets all reading frames for that exact shift
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

# INPUT: sequence
# OUTPUT: frames - all possible reading frames (list)
# get_all_possible_reading_frames creates a list from small lists of frames, for 3 shifts.
def get_all_possible_reading_frames(sequence):
    tmp_frames_list=[]
    frames=[]
    for shift in range(3):
        tmp_frames_list = reading_frame_shift(shift, sequence)
        frames = frames + tmp_frames_list

    return frames

# INPUT: frame_list - list of frame lists
# OUTPUT: new_list - list of frame lists (without lists with value None)
# remove_none_from_list function is needed because reading_frame function 
# creates None valued frames (if forming a frame was not possible)
def remove_none_from_list(frame_list):
    new_list = list(filter(lambda x: x[0] is not None, frame_list))
    return new_list

# INPUT: frames - reading frames list, strand - "+" or "-"
# sequence, total_genome_len - length of the genome
# OUTPUT: possible_orfs (list of dictionaries)
# find_possible_orfs - creates a dictionary item of found possible ORF
# It checks if start codon is right, does it have a stop codon.
# Calculates genomic position of the ORF
# Also, I've decided thats its better to find all possible candidates for ORFs.
# So after finding just one stop codon, code doesn't stop, and searches for the next stop codon.
# also, for "-" strand, it finds the possition not in its own sequence but in the starting genome
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

# INPUT: orfs_list (list), max_allowed_overlap - how much of the overlap is accepted
# OUTPUT: clean_orfs (list) of orfs with removed orfs which overlapped, sorted by start position
# filter_overlapping_orfs checks if there are overlapping ORFs, with a value which is accepted
# If overlapping is more than 15, it removes the ORF
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

# INPUT: orfs_lsit (list)
# OUTPUT: correct_orfs (list) - list of orfs which checks out for the length requirement
# filter_wrong_length_orfs - checks if the possible orf is the right length
# https://pmc.ncbi.nlm.nih.gov/articles/PMC3454372/
def filter_wrong_length_orfs(orfs_list):
    correct_orfs = []
    for orf in orfs_list:
        if 100 <= orf["length"] <= 700:
            correct_orfs.append(orf)
            
    return correct_orfs

# INPUT: orfs_list (list)
# OUTPUT: correct_orfs (list) - list of orfs i deemed worthy
# filtering_orfs joins two functions to check if possible orf is really an ORF
def filtering_orfs(orfs_list):
    correct_orfs = filter_wrong_length_orfs(orfs_list)
    correct_orfs = filter_overlapping_orfs(correct_orfs, max_allowed_overlap=15)

    return correct_orfs


# INPUT: STDIN
# OUTPUT: STDOUT
# main function joins every function and also using argparse creates a help page with choosable options
# for different functionalities of the program
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