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
    #for sequence in SeqIO.parse("input/JCVI-syn3.0.fasta", "fasta"):
    for sequence in SeqIO.parse("input/test.fasta", "fasta"):
        # print(sequence.id)
        # print(sequence.seq)
        # print(repr(sequence.seq))
        # print(len(sequence))
        seq = sequence.seq

    return str(seq)

def reverse_sequence(sequence):
    reverse_seq = []
    for n in list(sequence):
        if n == "A": reverse_seq.append("T")
        if n == "T": reverse_seq.append("A")
        if n == "C": reverse_seq.append("G")
        if n == "G": reverse_seq.append("C")
        # This doesn't work:
        #else: reverse_seq.append("C")

    reverse_seq.reverse()
    return "".join(reverse_seq)

def reading_frame(start, sequence):
    # Maybe need a check for start to be a negative number

    if start+3 <= len(sequence) and start >= 0:
        frame = sequence[start:start+3]
    else:
        return None
    return frame

def reading_frame_shift(shift, sequence):
    # iteruoti per butent ta shift
    shifted_frame = []
    for i in range(shift, len(sequence), 3):
        #print(i)
        shifted_frame.append(reading_frame(i, sequence))

    return shifted_frame


def get_all_possible_reading_frames(sequence):
    frames=[]
    # reikia poslinkiu 0,1,2
    # reikia indekso keitimo kas 3

    for shift in range(3):
        frames.append(reading_frame_shift(shift, sequence))

    return frames

def remove_none_from_list(frame_list):
    # THis doesn't work:
    # for i in reversed(frame_list):
        # if i is None:
        #     del frame_list[frame_list.index(i)]

    # Maybe lambda
    new_list = list(filter(lambda x: x != None, frame_list))
    # Gal reiktu rasyti x is not None
    return new_list

def join_reverse_forward_reading_frames():
    input_seq = fasta_in()
    reverse_seq = reverse_sequence(input_seq)

    forward = get_all_possible_reading_frames(input_seq)
    reverse = get_all_possible_reading_frames(reverse_seq)

    all_reading_frames = []
    for f, r in zip(forward, reverse):
        # Pasirodo negaliu naudoti append, nes prideda kaip sarasa
        all_reading_frames.extend(f)
        all_reading_frames.extend(r)
    #all_reading_frames = forward + reverse
    clean_frames = remove_none_from_list(all_reading_frames)
    return clean_frames


if __name__ == '__main__':
    # Original sequence
    #input_seq = fasta_in()
    #reverse_seq = reverse_sequence(input_seq)

    # Test for complementarity:
    #print(input_seq[0:10])
    #print(reverse_seq[0:10])

    # Tests for reading_frame function
    # print(input_seq[0:3])
    # print(reading_frame(0, input_seq))

    # Test for get_all_possible_reading_frames function
    # frames = get_all_possible_reading_frames(input_seq)
    # for i in range(len(input_seq)):
    #     print(frames[i])

    #forward_frames = get_all_possible_reading_frames(input_seq)
    #for frame in forward_frames: print(frame)

    orfs = join_reverse_forward_reading_frames()
    print(orfs)


    #ATGAAATAG