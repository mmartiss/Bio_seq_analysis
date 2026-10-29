from Bio import SeqIO

'''
How I visualize the program:

read_fasta()

define_states()

initialize_hmm()

train_hmm()

viterbi()

extract_genes()

evaluate_model()


'''

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

if __name__ == '__main__':
    seq = fasta_in("input/test_2.fasta", True)

    print(repr(seq))