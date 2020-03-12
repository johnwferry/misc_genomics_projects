from compsci260lib import *
from generate_HMM import generate_HMM


def summarize_sequence(state_sequence):
    """
    Given a state sequence, report the average length in each state

    Arguments:
        state_sequence (list of strings): generated sequence of states
    """

    # YOUR CODE GOES HERE


def main():
    """Load the sequence HMM file, then report the average length in each
    state"""
    seq_length = 1000
    state_sequence, observation_sequence = generate_HMM("HMM.sequence.txt",
                                                     seq_length)
    summarize_sequence(state_sequence)


if __name__ == '__main__':
    """Main method call, do not modify"""
    main()
