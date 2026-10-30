# HMM:

PDF is in the same directory.


# ORFs python program:

PDF is in the same directory.

## USED CRITERIA FOR POSSIBLE ORFs
Length : from 100 to 700 of codons \
Overlapping: max 15 \
Start codons: ATG, GTG, TTG \
Stop codons: TAA, TGA, TAG

For this genome: \
https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_037044595.1/ \
**Mycoplasma mycoides JCVI-syn3.0** \
NCBI ORF_finder results: 5879 ORFs \
My ORFs: 429 ORFs

## My choices for the program
### Step by step filtering
1. Search for START and STOP codons. If you have enough ram, you can search for all STOP codons for one start codon (just to see, that maybe the ORF doesn't end on the first seen STOP codon). Check if START and STOP codon is in the same reading frame.
2. Filter out by length.
3. Filtering by overlapping. Sort them by length. If the candidate does not exceed 15nt overlap, it is left alone, if not, filtered out.


### Length of the ORFs
I chose possible ORFs to be from 100 to 700 codon length, because of this article: \
https://pmc.ncbi.nlm.nih.gov/articles/PMC3454372/

I am not sure if this is correct way of doing it, but I believe that this is correct.

### Overlapping
I selected a maximum overlap of 15 nucleotides as a heuristic to reduce overlapping predictions while allowing limited overlap. Longer candidates are given priority. This rule may exclude genuine overlapping genes.

### Start and Stop codons
Start: ATG, GTG, TTG \
Stop: TAA, TGA, TAG

### Programming decisions
1. Used SeqIO for Fasta input.
2. Search for reading frames is a bit tedious, but it works. 
3. I calculate reverse strand position in the primary sequence. Like the genomic start and end.
4. Created a list of dictionary for every ORF, more information better? But it need a lot of RAM (not that much).
5. Argparse is not needed, but it helps using the program.
6. Most important, I have the debug mode, so maybe it helps to understand how the program works or what it does step by step.
7. I hope comments in the code makes everything clear.


### Disclaimer
Didn't remember how to use SeqIO. So kinda copied it from: \
https://lanadominkovic.medium.com/bioinformatics-101-reading-fasta-files-using-biopython-501c390c6820