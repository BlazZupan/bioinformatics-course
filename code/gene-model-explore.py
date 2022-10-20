from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt

organism_id = {
    "Hs_21q": "BA000005.3",  # Homo sapiens genomic DNA, chromosome 21q
    "Pl": "NC_001416.1",  # Phage lambda
    "Mt": "AL123456.3",  # Mycobacterium tuberculosis
    "Mg": "NC_000908.2",  # Mycoplasma genitalium
    "Ec": "NC_000913",  # E coli
}
organism = "Mg"

Entrez.email = "blaz.zupan@fri.uni-lj.si"
with Entrez.efetch(
    db="nucleotide",
    rettype="gbwithparts",
    retmode="text",
    id=organism_id[organism]
) as handle:
    rec = SeqIO.read(handle, "gb")

print("%s with %i features" % (rec.id, len(rec.features)))
print(rec.id)
print(rec.seq[:10])
print(len(rec.features))

cds = [f for f in rec.features if f.type == "CDS"]

# cds[0]
s = str(rec.seq)
s[685:688] <- ATG

# cds[-1]
r = str(rec.reverse_complement().seq)
s[580030:580033] <- CAT --rc--> ATG
r[len(s)-580033:len(s)-580033+3]

# len(cds) for c in cds

plt.cla()
plt.hist([len(c)/3 for c in cds], 50)
plt.xlabel("cds length [codons]")
plt.ylabel("frequency")
plt.title("Distribution of CDS lengths")
plt.savefig("cds.pdf")  # use pdf for inclusion in the report
plt.close()