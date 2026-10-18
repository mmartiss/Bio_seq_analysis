from Bio import SeqIO


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
def fasta_in():
    for sequence in SeqIO.parse("input/JCVI-syn3.0.fasta", "fasta"):
        # print(sequence.id)
        # print(sequence.seq)
        # print(repr(sequence.seq))
        # print(len(sequence))
        seq = sequence.seq

    return seq

if __name__ == '__main__':
    # Original sequence
    input_seq = fasta_in()
    print(input_seq)

