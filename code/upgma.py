from itertools import combinations, product
from collections.abc import Iterable


def average(xs):
    return sum(xs)/len(xs)


def flatten(xs):
    for x in xs:
        if isinstance(x, Iterable) and \
                not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x


class UPGMA:
    def __init__(self, file_name):
        """Read the distance matrix, initialize the clustering"""
        self.d = {}  # pairwise distances
        self.joins = []  # cluster joins in the order of joining
        self.clusters = None  # resulting clustering structure
        with open(file_name, "rt") as f:
            s = f.readline()
            self.names = s.strip().split("\t")
            for s in f.readlines():
                dis = s.strip().split("\t")
                self.d.update({(dis[0], n): float(d)
                               for n, d in zip(self.names, dis[1:])
                               if dis[0] != n})

    def cluster_distance(self, c1, c2):
        """
        Compute distance between two clusters.
        Example call: self.cluster_distance(("A", "B", "C"), ("D", "E"))
        """
        s = sum(self.d[a, b] for a, b in product(c1, c2))
        return s / (len(c1) * len(c2))

    def closest_clusters(self, clusters):
        """Return the pair of the closest clusters and their distance."""
        dis, pair = min((self.cluster_distance(list(flatten(c1)),
                                               list(flatten(c2))),
                         (c1, c2))
                        for c1, c2 in combinations(clusters, 2))
        return dis, pair

    def run(self):
        """Agglomerative hierarchical clustering"""
        clusters = [(name,) for name in self.names]
        while len(clusters) > 1:
            dis, (ca, cb) = self.closest_clusters(clusters)
            self.joins.append((ca, cb, dis))
            clusters = [c for c in clusters
                        if c not in [ca, cb]] + [(ca, cb)]
        self.clusters = clusters


if __name__ == "__main__":
    upgma = UPGMA("upgma-wiki-distances.txt")
    upgma.run()
    print(upgma.clusters)
