from collections import Counter
from bwt import *
from compsci260lib import *


def find(query, bwt_data):
    """Given a query sequence and a series of data structures containing
    various information about the reference genome, return a list containing
    all the locations of the query sequence in the reference genome.

    Args:
        query (str): query sequence to identify in the reference genome

    Returns:
        (list of ints): of all locations of the query sequence in the
                        reference genome
    """

    bwt, suffix_array, ranks, counts = bwt_data
    length = len(bwt)
    results = []
    query_value = query[len(query)-1]
    count_start_index = counts.get(query_value)
    count_end_index = length
    if query_value == 'A':
        next_value = 'C'
        count_end_index = counts.get(next_value)
    if query_value == 'C':
        next_value = 'G'
        count_end_index = counts.get(next_value)
    if query_value == 'G':
        next_value = 'T'
        count_end_index = counts.get(next_value)
    a = count_start_index + 1
    b = count_end_index

    qi = len(query) - 2
    while qi >= 0:
        qv = query[qi]
        a = counts.get(qv) + ranks.get(qv)[a - 1]
        b = counts.get(qv) + ranks.get(qv)[b-1]
        if a > b:
            return results
        qi -= 1
    """
    After the fnal step of the iteration (in which we have nished
considering the largest possible query sufx, i.e. the entire query sequence itself), if we nd that the range
of sorted suxes contains at least one row with a match to the query, we can use the range indices to obtain
the actual locations of those matches in the reference genome, locations which are themselves stored in the
sufx array
    
    """

    temp = []
    for i in range(a, b+1):
        if bwt[i] == query[qi+1]:
            temp.append(suffix_array[i])
    results = temp

    return sorted(results)

# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.
def rank(bwt_seq):
    """Takes as input a string transformed by the BWT. Returns a dictionary
    with characters as keys and lists as values. Each list contains the total
    number of occurrences for the corresponding character up until each
    position in the BWT-transformed string (i.e., its rank).

    For example:
        rank('ACTGA$TA')['$'] --> [0, 0, 0, 0, 0, 1, 1, 1]
        rank('ACTGA$TA')['A'] --> [1, 1, 1, 1, 2, 2, 2, 3]
        rank('ACTGA$TA')['C'] --> [0, 1, 1, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['G'] --> [0, 0, 0, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['T'] --> [0, 0, 1, 1, 1, 1, 2, 2]

    Args:
        bwt_seq (str): BWT-transformed sequence

    Returns:
        (dict): with characters as keys and lists of integers
        containing the total number of occurrences for the
        corresponding character up until each position in the
        BWT-transformed string (i.e., its rank)
    """
    rank = {}
    characters = set(bwt_seq)
    for character in characters:
        rank[character] = [0]
    rank[bwt_seq[0]] = [1]
    for letter in bwt_seq[1:]:
        for k, v in list(rank.items()):
            v.append(v[-1] + (k == letter))
    return rank


def make_suffix_array(seq):
    """Makes the suffix array of a given input string sequence.

    For example:
        make_suffix_array('GATTACA$') --> [7, 6, 4, 1, 5, 0, 3, 2]

    Args:
        seq (str): input string with an EOF character

    Returns:
        (list): of integers of the suffix array of the input string.
    """
    suffixes = {}
    for x in range(len(seq)):
        suffixes[seq[x:]] = x
    suffix_array = [suffixes[suffix] for suffix in sorted(suffixes.keys())]
    return suffix_array


def count_smaller_chars(seq):
    """Takes as input a string. Returns a dictionary with characters as keys
    and integers as values. The integers track the number of characters in the
    input string which are lexicographically smaller than the corresponding
    character key.

    For example, using an input DNA sequence like 'GATTACA':
        count_smaller_chars('GATTACA')['A'] --> 0
            (A, being lexicographically first in a DNA sequence,
            should always return 0)

        count_smaller_chars('GATTACA')['C'] --> 3
            (C, being second, should return the number of A's, which here is 3)

        count_smaller_chars('GATTACA')['G'] --> 4
            (G, being third, should return the number of A's or C's,
            which here is 4)

        count_smaller_chars('GATTACA')['T'] --> 5
            (T, being fourth, should return the number of A's or C's or G's,
            which here is 5)
    """
    characters = set(seq)
    cntr = Counter(seq)
    total = 0
    counts = {}
    for character in sorted(characters):
        counts[character] = total
        total += cntr[character]
    return counts


def make_all(reference):
    """Takes as input a reference string. Returns the data structures necessary
    to perform efficient exact string matching searches.

    Args:
        reference (str): reference string to create data structures for

    Returns:
        tuple of
        (str) forward bwt of the reference string
        (list of int): suffix_array of the reference string
        (list of int): ranks of the forward bwt
        (dict of str to int): smaller character counts
    """
    counts = count_smaller_chars(reference)
    reference = reference + '$'
    suffix_array = make_suffix_array(reference)
    bwt = forward_bwt(reference)
    ranks = rank(bwt)
    return bwt, suffix_array, ranks, counts


if __name__ == '__main__':
    # example query sequence
    query_sequence = "AAACGA"
    # example reference sequence
    sequence = "AAAAAAAAACGATAGAGA"
    find(query_sequence, make_all(sequence))

