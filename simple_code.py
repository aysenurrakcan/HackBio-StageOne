# Ayşenur Akcan
# HackBio Stage 1
# This script translates a DNA sequence to a protein
# and calculates the Hamming distance between two strings.

from Bio.Seq import Seq

def translate_dna(dna_seq):
    """
    Translates a given DNA sequence into a protein sequence
    using the Standard codon table.
    Includes error handling for sequences not divisible by 3.
    """
    dna_seq = dna_seq.upper()
    
    # Check sequence length validity
    if len(dna_seq) % 3 != 0:
        raise ValueError("DNA sequence length must be divisible by 3 for proper translation.")
    
    dna = Seq(dna_seq)
    # Specify codon table explicitly
    protein = dna.translate(table="Standard")
    return protein


def hamming_distance(a: str, b: str, pad_char: str = "_") -> int:
    """
    Calculates the Hamming distance between two strings.
    Pads shorter string with a character if lengths differ.
    """
    a = a.lower()
    b = b.lower()

    max_len = max(len(a), len(b))
    a_padded = a.ljust(max_len, pad_char)
    b_padded = b.ljust(max_len, pad_char)

    return sum(ch1 != ch2 for ch1, ch2 in zip(a_padded, b_padded))


# --- Example usage ---

# DNA to protein translation
dna_sequence = "ATCTGGCCATAATTCGTTGGCTGAAAGGAGTGCGCCTCCGATAG"
try:
    protein = translate_dna(dna_sequence)
    print("Protein sequence:", protein)
except ValueError as e:
    print("Error:", e)

# Hamming distance calculation
slack = "Ayşenur Akcan"
twitter = "aysenura"

distance = hamming_distance(slack, twitter)
print("Hamming distance:", distance)
