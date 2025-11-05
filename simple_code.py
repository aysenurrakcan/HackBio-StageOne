%pip install biopython
from Bio.Seq import Seq

dna = Seq("ATCTGGCCATAATTCGTTGGCTGAAAGGAGTGCGCCTCCGATAG")
protein = dna.translate()
print(protein)

output:
IWP*FVG*KECASD


#HAMMING DISTANCE
def hamming_distance(a: str, b: str, pad_char: str = "_") -> int:
   
    a = a.lower()
    b = b.lower()

    max_len = max(len(a), len(b))
    a_padded = a.ljust(max_len, pad_char)
    b_padded = b.ljust(max_len, pad_char)

    return sum(ch1 != ch2 for ch1, ch2 in zip(a_padded, b_padded))



slack = "Ayşenur Akcan"
twitter = "aysenura"

distance = hamming_distance(slack, twitter)
print("Hamming distance:", distance)
