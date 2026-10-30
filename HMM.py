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

    states = ["ig", "start1", "start2", "start3", "e1", "e2", "e3", "i0", "i1", "i2", "s1", "s2", "s3"]

    transitions = {
        "ig": {"ig": 0.99, "start1": 0.01},

        "start1": {"start2": 1.00},
        "start2": {"start3": 1.00},
        "start3": {"e1": 1.00},

        "e1": {"e2": 0.95, "i1": 0.05},
        "e2": {"e3": 0.95, "i2": 0.05},
        "e3": {"e1": 0.84, "i0": 0.05, "s1": 0.11},

        "i0": {"i0": 0.99, "e1": 0.01},
        "i1": {"i1": 0.99, "e2": 0.01},
        "i2": {"i2": 0.99, "e3": 0.01},

        "s1": {"s2": 1.00},
        "s2": {"s3": 1.00},
        "s3": {"ig": 1.00}
    }

    emissions = {
        "ig": {"A": 0.30, "C": 0.20, "G": 0.20, "T": 0.30},

        "start1": {"A": 1.0},
        "start2": {"T": 1.0},
        "start3": {"G": 1.0}
    }

    for state in states:
        total = sum(transitions[state].values())
        if not abs(total - 1.0) < 1e-9:
            print(f"Invalid transitions for {state}: {total}")

    print(repr(seq))


if __name__ == '__main__':
    main()