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
        # I WANT TO RETURN NOTHING, THIS BRAKES CODEEE
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
    # THis doesn't work:
    # for i in reversed(frame_list):
        # if i is None:
        #     del frame_list[frame_list.index(i)]

    # Maybe lambda
    #new_list = list(filter(lambda x: x != None, frame_list))
    # Gal reiktu rasyti x is not None

    #Nebeveikai, kai padariau sarasiukais
    new_list = list(filter(lambda x: x[0] is not None, frame_list))
    return new_list

# def join_reverse_forward_reading_frames():
#     input_seq = fasta_in()
#     reverse_seq = reverse_sequence(input_seq)

#     forward = get_all_possible_reading_frames(input_seq)
#     reverse = get_all_possible_reading_frames(reverse_seq)

#     # all_reading_frames = []
#     # for f, r in zip(forward, reverse):
#     #     # Pasirodo negaliu naudoti append, nes prideda kaip sarasa
#     #     all_reading_frames.extend(f)
#     #     all_reading_frames.extend(r)

#     # for f in forward:
#     #     all_reading_frames.extend(f)
#     # for r in reverse:
#     #     all_reading_frames.extend(r)


#     #all_reading_frames = forward + reverse
#     #clean_frames = remove_none_from_list(all_reading_frames)
#     return forward, reverse

def find_possible_orfs(frames, strand, sequence):
    which = 0
    possible_orfs = []
    for frame in frames:
        if frame[0] == 'ATG' or frame[0] == 'GTG' or frame[0] == 'TTG':
            for stop in frames:
                if frame[1] == stop[1] and frame[2] < stop[2] and (stop[0] == 'TAA' or stop[0] == 'TAG' or stop[0] == 'TGA'):
                    seq = sequence[frame[2]:stop[2]+3]
                    length = len(seq)//3 # nes antraip float
                    orf = {
                        "start": frame[0],
                        "stop": stop[0],
                        "shift": frame[1],
                        "strand": strand,
                        "sequence": seq,
                        "length": length,
                        "index" : which
                    }
                    possible_orfs.append(orf)
                    which += 1

    return possible_orfs



if __name__ == '__main__':
    # Original sequence
    input_seq = fasta_in()
    reverse_seq = reverse_sequence(input_seq)
    forward = get_all_possible_reading_frames(input_seq)
    reverse = get_all_possible_reading_frames(reverse_seq)
    forward_clean_frames = remove_none_from_list(forward)
    reverse_clean_frames = remove_none_from_list(reverse)

    print(forward_clean_frames)
    print(reverse_clean_frames)

    f_orfs = find_possible_orfs(forward_clean_frames, "+", input_seq)
    r_orfs = find_possible_orfs(reverse_clean_frames, "-", reverse_seq)


    print(f_orfs)
    print(r_orfs)
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

    # frames = join_reverse_forward_reading_frames()
    # print(frames)

    # forward, reverse = join_reverse_forward_reading_frames()

    # print(forward)
    # print(len(frames))
    #print(len(input_seq))

    #ATGAAATAG