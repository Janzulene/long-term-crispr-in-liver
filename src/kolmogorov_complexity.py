from typing import Generator
from itertools import product
from collections import defaultdict

import numpy as np


def get_kmer_list(k: int) -> Generator[str, None, None]:
    "Generate all possible kmers of length k."
    return ("".join(p) for p in product("ATCG", repeat=k))


def calc_kmer_entropy(seq: str, k: int) -> float:
    """Calculate the normalized entropy of kmers in a sequence.

    Args:
        seq (str): DNA sequence.
        k (int): the length of kmer.

    Returns:
        float: the normalized entropy of kmers.
    """
    kmer_dist = {i: 0 for i in get_kmer_list(k)}

    for i in range(len(seq) - k + 1):
        # Skip non-ACGT bases
        if seq[i:i + k] in kmer_dist:
            kmer_dist[seq[i:i + k]] += 1

    count = np.array(list(kmer_dist.values()))
    prob  = count / count.sum()
    prob  = prob[prob > 0]
    return -np.sum(prob * np.log2(prob)) / np.log2(len(kmer_dist))


def calc_kolmogorov_complexity(seq: str, max_k: int = 4) -> float:
    """Calculate a proxy for the Kolmogorov complexity of a sequence.

    Defined here as the sum of normalized kmer entropies for k from 2 to max_k.

    Args:
        seq (str): DNA sequence.
        max_k (int, optional): the largest kmer length to consider, must be
            greater than 2. Defaults to 4.

    Returns:
        float: the sequence complexity score.
    """
    assert max_k > 2, "max_k should be greater than 2!"
    seq = seq.upper()
    return sum(calc_kmer_entropy(seq, k) for k in range(2, max_k + 1))
