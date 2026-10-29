from Bio import SeqIO

'''
How I visualize the program:

read_fasta() +

define_states()

initialize_hmm()

train_hmm()

viterbi()

extract_genes()

evaluate_model()
'''

# Returns only the last sequence in the fasta file
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



def main():
    seq = fasta_in("input/test_2.fasta", True)

    states = ["ig", "e1", "e2", "e3", "i"]

    transitions = {
        "ig":
        "e1":
        "e2":
        "e3":
        "i":
    }

    print(repr(seq))


if __name__ == '__main__':
    main()