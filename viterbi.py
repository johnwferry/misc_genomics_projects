
import sys
from math import log
from compsci260lib import *


def viterbi_decoding(input_file, hmm_file):
    """Calculate the viterbi decoding of an input sequence

    Arguments:
        input_file (str): path to input fasta file
        hmm_file (str): path to HMM file

    Returns:
        A list of dictionaries of segments in each state. An example output may
        look like:

        [
            {‘start’: 0, ‘end’: 12, ‘state’: ‘state2’},
            {‘start’: 13, ‘end’: 20, ‘state’: ‘state1’},
            ...
        ]
    """

    # open hmm file
    try:
        f_hmm_file = open(hmm_file, 'r')
    except IOError:
        print(("IOError: Unable to open HMM file: %s." % hmm_file))
        print("Exiting.")
        sys.exit(1)

    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)

    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]

    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row

    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()

    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row

    f_hmm_file.close()

    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1

    print("Done reading sequence of length ", len(emit_str))

    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)]

    # Initialize the first column of the matrix
    for i in range(K):
        in_index = get_emit_index(emit_str[0].upper(), emitted_symbols)
        viterbi[i][0] = log(emit_probs[i][in_index]) + log(initial_probs[i])

    # Build the matrix column by column
    for j in range(1, len(emit_str)):
        in_index = get_emit_index(emit_str[j].upper(), emitted_symbols)

        for i in range(K):
            # Compute the entries viterbi[i][j] and pointers[i][j]
            # Tip: Use float('-inf') for the value of negative infinity
            """YOUR CODE HERE"""
            maxv = float('-inf')
            maxi = -1
            for k in range(K):
                temp = log(transitions[k][i]) + viterbi[i][j-1]
                if temp > maxv:
                    maxv = temp
                    maxi = k
            viterbi[i][j] = log(emit_probs[i][in_index]) + maxv
            pointers[i][j] = maxi
    maxend = float('-inf')
    maxk = -1
    for i in range(K):
        temp = viterbi[i][len(emit_str)-1]
        if temp > maxend:
            maxend = temp
            maxk = i





    # Traceback, stored as a list segment lengths in each state as dictionaries
    """YOUR CODE HERE"""
    traceback = []
    j = len(emit_str)-1
    endpointer = pointers[maxk][j]
    end = j
    seglength = 1
    while j > 0:
        j = j - 1
        temp = pointers[endpointer][j]
        if temp is not endpointer:
            traceback.append({'start': j, 'end': end, 'state': 'state'+str(endpointer)})
            endpointer = temp
            end = j
            seglength = 1
        else:
            seglength += 1
    traceback.append({'start': j, 'end': end, 'state': 'state'+str(endpointer)})

    return traceback


def count_segments(vit_ret):
    """Calculate the number of segments appearing in each state of
    the viterbi path

    Arguments:
        vit_ret (list of dicts): dictionary of lengths in each state.
            see: return value of viterbi_decoding

    Returns:
        a dictionary of states to number of occurrences in the state. e.g.
        {'state0': 10, 'state1': 2}
    """

    """YOUR CODE HERE"""
    dictret = {}
    for l in vit_ret:
        if l.get('state') not in dictret:
            dictret[l.get('state')] = 0
        dictret[l.get('state')] = dictret[l.get('state')] + 1
    return dictret

def get_emit_index(input_val, alphabet):
    """Get the index of the emission value in the alphabet

    Note: This will be a useful function for indexing the emission
    probabilities"""
    return alphabet.index(input_val)


def main():
    hmm_file = 'HMM.methanococcus.txt'
    input_file = 'artificial.genome.fasta'
    vit_ret = viterbi_decoding(input_file, hmm_file)

    # report the number of segments that exist in each state.
    counts = count_segments(vit_ret)
    """YOUR CODE HERE"""


if __name__ == '__main__':
    """Call main(), do not modify."""
    main()
