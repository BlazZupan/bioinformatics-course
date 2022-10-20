import os.path
from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "blaz.zupan@fri.uni-lj.si"

organism_id = {
    "Hs_21q": "BA000005.3",  # Homo sapiens genomic DNA, chromosome 21q
    "Pl": "NC_001416.1",  # Phage lambda
    "Mt": "AL123456.3",  # Mycobacterium tuberculosis
    "Mg": "NC_000908.2",  # Mycoplasma genitalium
    "Ec": "NC_000913",  # E coli
}
organism = "Mg"

filename = f"data/{organism}.fasta"
if os.path.exists(filename):
    print(f"Loading {organism} sequence from local file.")
    data = SeqIO.read(filename, "fasta")
else:
    print(f"Fetching {organism} ...")
    with Entrez.efetch(
        db="nucleotide",
        id=organism_id[organism],
        rettype="gbwithparts",
        retmode="text"
    ) as handle:
        data = SeqIO.read(handle, format="genbank")

    print("Storing ...")
    with open("data/%s.fasta" % organism, "w") as f:
         SeqIO.write([data], f, "fasta")

freq = Counter(str(data.seq))
prob = {n: freq[n]/len(data) for n in freq}
