%pip install biopython
from Bio.Seq import Seq

dna = Seq("ATCTGGCCATAATTCGTTGGCTGAAAGGAGTGCGCCTCCGATAG")
protein = dna.translate()
print(protein)

output:
IWP*FVG*KECASD
