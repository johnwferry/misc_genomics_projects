from compsci260lib import *
import sys
import random
import os


def generate_HMM(hmm_file, seq_length):
    """ Load an HMM specification from file, and generate state and observed
        sequences through sampling of the HMM.

        The HMM will have states labeled as strings and observations as (string)
        characters which will be useful when we generate an HMM with nucleotide
        sequence observations:

        An example return sequence for an occasionally dishonest casino HMM
        with dice rolls may look like:
            (['F', 'L', 'L'], ['2', '6', '6'])

        Arguments:
            hmm_file (str): path to the HMM file
            seq_length (int): the length of the sequence we will generate

        Returns:
            a tuple of
            (the state sequence as strings,
             observed sequence as single character strings)
    """
    if not os.path.exists(hmm_file):
        raise ValueError("Can't open HMM parameter file: %s" % hmm_file)

    f = open(hmm_file, "r")

    # read the state names
    states = f.readline().strip().split()

    # read the initial probabilities
    initial_probs = f.readline().strip().split()
    initial_probs = [float(p) for p in initial_probs]

    # read and store the transition matrix
    # NOTE: this is stored as a dictionary of lists, so the rows of the
    #       matrix are "named" with the state name as a key, while the
    #       columns of the matrix correspond to indices into the list
    transitions = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        transitions[state] = [float(p) for p in matrix_row]

    # read in all the observable characters
    observable_chars = f.readline().strip().split()

    # read the emission matrix
    # NOTE: this is stored as a dictionary of lists, so the rows of the
    #       matrix are "named" with the state name as a key, while the
    #       columns of the matrix correspond to indices into the list
    emission = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        emission[state] = [float(p) for p in matrix_row]
        # normalize
        sum_emissions = sum(emission[state])
        emission[state] = [x / sum_emissions for x in emission[state]]
    f.close()

    ran = random.random()
    def get_random_result(li):
        temp = random.random()
        runningtotal = 0
        for x in range(0,len(li)):
            runningtotal += li[x]
            if temp < runningtotal:
                return x
        return x
    def to_nuc(input):
        if input == 0:
            return 'A'
        if input == 1:
            return 'C'
        if input == 2:
            return 'G'
        return 'T'


    first_state = get_random_result(initial_probs)
    sl = []
    sl.append(first_state)
    el = []
    el.append(to_nuc(get_random_result(emission[states[sl[0]]])))
    for i in range(1, seq_length):
        sl.append(get_random_result(transitions[states[sl[i-1]]]))
        el.append(to_nuc(get_random_result(emission[states[sl[i]]])))

    for i in range(0, seq_length):
        sl[i] = sl[i] + 1

    s1posc = []
    s2posc = []
    cnt = 1
    pp = sl[0]
    for i in range(1, seq_length):
        if sl[i] == pp:
            cnt += 1
        else:
            if pp == 1:
                s1posc.append(cnt)
                cnt = 1
            if pp == 2:
                s2posc.append(cnt)
                cnt = 1
            pp = sl[i]

    for i in range(0, seq_length):
        print(sl[i], end='')
    print('')
    for i in range(0, seq_length):
        print(el[i], end='')
    print('')

    print("Avg state 1 length: " + str(sum(s1posc)/len(s1posc)))
    print("Avg state 2 length: " + str(sum(s2posc)/len(s2posc)))
    return (sl, el)  # a tuple containing the state and observation sequences


def audit_casino(state_sequence, observed_sequence, subsequence):
    """Given state and observed sequences, report the state sequences and
    frequency of those state sequences that emitted the rolls represented
    in subsequence  
    """
    dic = {}
    l = len(subsequence)
    for i in range(0,len(state_sequence)-l):
        fail = False
        for j in range(0,l):
            if subsequence[j] != observed_sequence[j+i]:
                fail = True
        if not fail:
            s = ''
            for k in state_sequence[i:i+l]:
                s = s + k
            if s in dic:
                dic[s] = dic[s] + 1
            else:
                dic[s] = 1
    print(dic)



def main():
    """Generate state and observed sequences and report the frequency of state
    sequences emitting 1662"""
    seq_length = 100000
    state_sequence, observed_sequence = generate_HMM("HMM.sequence.txt",
                                                     seq_length)
    audit_casino(state_sequence, observed_sequence, ['1','6','6','2'])


if __name__ == '__main__':
    """Main method call, do not modify"""
    main()
