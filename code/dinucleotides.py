import random
from Bio import SeqIO
from collections import Counter
from operator import mul
from functools import reduce
from scipy.stats import binom
import numpy as np
from matplotlib import pyplot as plt


def product(lst):
    return reduce(mul, lst, 1)


def tuple_walk(s, k=2):
    """Generate subsequences of length k."""
    for i in range(len(s)-1):
        yield s[i:i+k]


def count_tuple(s, w):
    """Count words w in sequence s."""
    return sum(1 for ss in tuple_walk(s, len(w)) if ss == w)


organism = "Hs_21q"
# organism = "Mg"
filename = f"data/{organism}.fasta"
s = str(SeqIO.read(filename, "fasta").seq)

subsequence = "CG"
count = count_tuple(s, subsequence)
print(f"Count of {subsequence} on original sequence: {count}")

n = len(s)
p_nucleotide = {k: v/n for k, v in Counter(s).items()}
p = product(p_nucleotide[x] for x in subsequence)
print("Expected: %d" % (p * n))


def permutation_counts(seq, t, epochs=100):
    seq = list(seq)
    counts = []
    for i in range(epochs):
        print(i)
        random.shuffle(seq)
        counts.append(count_tuple("".join(seq), t))
    return counts

distribution = permutation_counts(s, subsequence)
p_value = sum(1. for x in distribution if x >= count)/len(distribution)
print("p-value: %5.4f" % p_value)

plt.cla()
plt.hist(distribution, 20)
plt.title(f"Permutation test, frequency of {subsequence}")
plt.xlabel(f"occurrences of {subsequence} in permuted genome")
plt.ylabel("frequency")
plt.savefig("res0.pdf")

# the right way to assess the distribution
# show also [(k, binom.pmf(k, 10, 0.3)) for k in range(10)]

n = len(s)
p_nucleotide = {k: v/float(n) for k, v in Counter(s).items()}
p = product(p_nucleotide[x] for x in cg)

# Percent point function (inverse of cdf - percentiles)
n = n - len(cg) + 1
x = np.arange(binom.ppf(0.01, n, p), binom.ppf(0.99, n, p))
binom.pmf(x, n-1, p)

# Probability mass function.
plt.plot(x, binom.pmf(x, n, p), 'bo', ms=1, label='binom pmf')
plt.savefig("res.pdf")
print("for results: open res.pdf")
plt.close()

# Cumulative density function.
print("p-value (from binomial distribution): %5.4f" %
      (1. - binom.cdf(count, n, p)))
