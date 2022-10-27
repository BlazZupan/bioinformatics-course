from collections import defaultdict


def score(a, b):
    """A scoring function."""
    return 1 if a == b else 0


def longest_common_substring(s, t):
    """Returns a dynamic programming table for the longest
    common substring of sequences s and t, and pointers to
    optimal previous predecessor cell for each of the cells
    in the table."""
    M = defaultdict(int)
    P = {}
    for i, si in enumerate(s):
        for j, tj in enumerate(t):
            M[i, j], P[i, j] = max(
                (M[i-1, j], (i-1, j)),
                (M[i, j-1], (i, j-1)),
                (M[i-1, j-1] + score(si, tj), (i-1, j-1))
            )
    return M, P


def trace_back(s, t, M, P):
    """Trace-back through the dynamic programming table to
    return the longest common subsequence."""
    i = len(s)-1
    j = len(t)-1
    marked_nodes = set()
    result = ""
    while M[i, j] != 0:
        if P[i, j] == (i - 1, j - 1):
            result = s[i] + result
            marked_nodes.add((i, j))
        i, j = P[i, j]
    return result, marked_nodes


def print_dp_table(M, s, t, marks=None):
    """Print dynamic programming table"""
    print("  " + " ".join(t))
    for i in range(len(s)):
        if marks:
            print("%-2s" % s[i] +
                  " ".join("%d" % M[i, j]
                           if (i, j) in marks
                           else " " for j in range(len(t))))
        else:
            print("%-2s" % s[i] +
                  " ".join("%d" % M[i, j] for j in range(len(t))))


s0 = "AACCTTGG"
s1 = "ACACTGTGA"

table, previous = longest_common_substring(s0, s1)
print(f"s = {s0}")
print(f"t = {s1}")
print("Longest common substring length:", table[len(s0) - 1, len(s1) - 1])

substring, marked_nodes = trace_back(s0, s1, table, previous)
print(f"Longest common substring: {substring}")
print()

print("Dynamic programming table:")
print_dp_table(table, s0, s1, marked_nodes)
