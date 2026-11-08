from Bio import SeqIO

import argparse
import sys

def fasta_in(file, debug):
    for sequence in SeqIO.parse(file, "fasta"):
        if debug:
            print(sequence.id)
            print(sequence.seq)
            print(repr(sequence.seq))
            print(len(sequence))
        seq = sequence.seq

    return str(seq)

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
        "-o", "--output", type=str, default="output.fasta", help="Output file name or path"
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

if __name__ == '__main__':
    main()