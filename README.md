# BioinformaticsTools
Collection of various helper methods for working with genomic data.

- ```count_nucleotides```:
    - Count A, C, G, and T or U bases in provided DNA/RNA sequence
- ```transcribe_seq```:
    - Convert provided DNA sequence to RNA sequence
- ```get_reverse_complement```:
    - Return reverse complement of provided DNA sequence
- ```get_GC_content```:
    - Calculate percentage of G and C nucleotides in provided DNA sequence
- ```read_fasta_seqs```:
    - Return dictionary of labels and sequences read from provided FASTA file
- ```get_highest_GC_content```:
    - Return label and value of sequence with highest GC content from provided data dictionary
- ```find_motif_occurences```:
    - Find all non-overlapping occurences of provided motif (pattern) in sequence
- ```one_hot_encode```:
    - One-hot encode provided DNA sequence as np.ndarray