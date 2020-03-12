from compsci260lib import *


def assembly_tester(reads_file, contig0_file, contig1_file):
    """
    Read in the contig files, verify the reads in each contig, and estimate
    the gap length between the contigs
    """

    # load the fasta files
    reads = get_fasta_dict(reads_file)
    contig0 = get_fasta_dict(contig0_file)
    contig1 = get_fasta_dict(contig1_file)

    # make sure the contigs pass the criteria in part d
    """YOUR CODE HERE"""

    # determine how the reads map to the contigs
    contig_reads = find_contig_reads(reads, contig0, contig1)

    # report the reads whose ends map to both contig0 and contig1
    # and their estimated gap lengths

    """YOUR CODE HERE"""


def find_contig_reads(reads, contig0, contig1):
    """
    Determine whether the sequencing reads map to contig0/contig1/both/neither
    and where in the contig it matches. `reads` will contain both ends 'a' and
    'b', but you will return a dictionary using the read name as a whole
    (without 'a' and 'b').

    It should return a dictionary mapping the name of the mated pair of reads to:
        - 'contig_a' (str): the contig in which read 'a' was found in as `contig0',
          `contig1', or None.
        - start_a (int): the start position (1-indexed) read end 'a' mapped to
          for its respective contig (None if not found in any contig)
        - end_a (int): the end position (1-indexed) read end 'a' mapped to
          for its respective contig (None if not found in any contig)
        - 'contig_b' (str): the contig in which read 'b' was found in as `contig0',
          `contig1', or None.
        - start_b (int): the start position (1-indexed) read end 'b' mapped to
          for its respective contig (None if not found in any contig)
        - end_b (int): the end position (1-indexed) read end 'b' mapped to
          for its respective contig (None if not found in any contig)

    The returned dictionary should look something like:
    {
        'seq1': {
            'contig_a': 'contig1',
            'start_a': 1,
            'end_a': 100,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        'seq2': {
            'contig_a': 'contig0',
            'start_a': 101,
            'end_a': 200,
            'contig_b': 'contig1'
            'start_b': 201,
            'end_b': 300,
        },
        'seq3' : {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        ...
    }

    Arguments:
        reads (dict str to str): dictionary of sequencing reads
        contig0 (dict str to str): dictionary of reads in contig0
        contig1 (dict str to str): dictionary of reads in contig1

        see: get_fasta_dict

    Returns:
        Dictionary mapping reads to information of their mapping contig.
    """

    """YOUR CODE HERE"""

    return {}


def main():
    """Call assembly tester with provided fasta files"""
    reads_file = 'paired.reads.fasta'
    contig0_file = 'contig0.fasta'
    contig1_file = 'contig1.fasta'
    assembly_tester(reads_file, contig0_file, contig1_file)


if __name__ == '__main__':
    """Call main(). Do not modify this block."""
    main()
