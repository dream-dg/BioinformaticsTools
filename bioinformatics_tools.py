import numpy as np

def count_nucleotides(input, rna=False):
    """
    Count A, C, G, and T or U bases in provided DNA/RNA sequence
    """
    if rna:
        t_u_count = input.count("U")
    else:
        t_u_count = input.count("T")
    
    return input.count("A"), input.count("C"), input.count("G"), t_u_count

def transcribe_seq(input):
    """
    Convert provided DNA sequence to RNA sequence
    """
    return input.replace("T", "U")

def get_reverse_complement(input):
    """
    Return reverse complement of provided DNA sequence
    """
    complements = {"A":"T", "C":"G", "G":"C", "T":"A"}
    rev_comp = ""
    for i in input[::-1]:
        rev_comp += complements.get(i)

    return rev_comp

def get_GC_content(input):
    """
    Calculate percentage of G and C nucleotides in provided DNA sequence
    """
    return (input.count("G") + input.count("C")) / len(input) * 100

def read_fasta_seqs(file_name):
    """
    Return dictionary of labels and sequences read from provided FASTA file
    """
    data = {}
    with open(file_name, "r") as f:
        for line in f:
            if line.startswith('>'):
                curr = line.strip()[1:]
                data[curr] = ""
            else:
                data[curr] += line.strip() 
    return data

def get_highest_GC_content(data):
    """
    Return label and value of sequence with highest GC content from provided data dictionary
    """
    highest_val = 0
    for key, seq in data.items():
        gc_content = get_GC_content(seq)
        if gc_content > highest_val:
            highest_key = key
            highest_val = gc_content

    return highest_key, highest_val

def find_motif_occurences(seq, motif):
    """
    Find all non-overlapping occurences of provided motif (pattern) in sequence
    """
    locations = []
    for i in range(len(seq) - len(motif) + 1):
        if motif == seq[i:i+len(motif)]:
            locations.append(str(i+1))
    return locations

def one_hot_encode(sequence, alphabet="ACGT", neutral="N", neutral_val=0.25, dtype=np.float32):
    """
    One-hot encode provided DNA sequence as np.ndarray
    """
    def to_uint8(string):
        return np.frombuffer(string.encode("ascii"), dtype=np.uint8)
    encoded_sequence = neutral_val*np.ones((len(sequence), len(alphabet)), dtype=dtype)
    hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
    hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
    hash_table[to_uint8(neutral)] = neutral_val
    hash_table = hash_table.astype(dtype)
    encoded_sequence[:len(sequence)] = hash_table[to_uint8(sequence.upper())]
    return encoded_sequence.T