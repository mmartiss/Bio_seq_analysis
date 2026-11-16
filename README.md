# Motif search

## 1. Motifs.py input/output

The program identifies potential DNA sequence motifs using two methods: **k-mer enrichment analysis** and **Gibbs sampling**.

**Input:** A FASTA file containing at least two DNA sequences consisting of A, T, C, and G nucleotides.
**Output:** Results are displayed in STDOUT and saved to a `.txt` file.

| Argument               | Description                         |
| ---------------------- | ----------------------------------- |
| `-h`, `--help`         | Display help message                |
| `-i`, `--input`        | Specify input FASTA file            |
| `-d`, `--default-file` | Use the default test file           |
| `--debug`              | Enable detailed debugging output    |
| `-o`, `--output`       | Specify output file                 |
| `-k`, `--kmer`         | Set motif length                    |
| `--seed`               | Set random seed                     |
| `--iter`               | Set Gibbs sampling iterations       |
| `--restarts`           | Set independent Gibbs sampling runs |

## 2. K-mer Analysis and Background Model

Here I extract all overlapping k-mers from the input seqs. Make sure not to miss any possible positions. Used `Counter` to count k-mer occurrences.

A background model is calculated from the overall frequencies of A, T, C, and G, assuming that nucleotides occur independently. These frequencies are used to estimate the expected number of occurrences for each k-mer.

The enrichment value is calculated as the ratio between observed and expected occurrences. 

I have sorted results by enrichment, not frequency. Maybe I should add an option to choose which sorting the person needs.
## 3. Gibbs Sampling Implementation

Select random motif positions in each seq. Used a loop to iterate (there is an option for how many times), select one random sequence to exclude. And then during each loop constructed position porbability profile for remaining motifs.

Added a pseudocount of `+1` for each nucleotide count, so now zero probabilities would be there. Then just for excluded seq evaluate all possible motif positions, using a ratio, between motif profile probabilities and background nucleotide probabilities.

Used `random.choices()` so algo could maybe see more alternative solutions, so not always locally preferred position is chosen (highest score)

The program evaluates solutions using the number of mismatches relative to the consensus sequence. Also added independent restarts, used with `--restarts` flag. 
## 4. AI disclamer

AI assistance was used primarily for output formatting, refactoring part of the Gibbs sampling loop into a separate function, and minor implementation suggestions.
All parts where AI was used are commented in the code, with what exactly was done.
Model used: GPT 6-sol light.


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
