from bwt import *
from read_aligner import *
from compsci260lib import *


def reverse_complement(seq):
    """Returns the reverse complement of the input string."""
    comp_bases = {'A': 'T',
                  'C': 'G',
                  'G': 'C',
                  'T': 'A'}
    rev_seq = list(seq)
    rev_seq = rev_seq[::-1]
    rev_seq = [comp_bases[base] for base in rev_seq]
    return ''.join(rev_seq)


def read_mapper(patient, bacteria):
    """For a given patient and bacterial species, return a vector which counts
    the number of the patient's reads which map to each of the locations in the
    reference genome for the bacterial species.

    Args:
        patient (str): The patient name. Can be either "patient1", "patient2",
                       or "patient3". This will be useful for loading your
                       patient prevalence files from disk as described in
                       part c.

        bacteria (str): of the bacteria name as named in the reference_genomes
                        folder. e.g. "Bacteroides_ovatus", "Vibrio_cholerae",
                        etc. This will be useful for loading the appropriate
                        reference genome fasta file from disk.

    Returns:
        (list of ints): vector of aligned read counts to the genome.
        i.e. [c_1, c_2, ..., c_i, ..., c_n], where n=length of the genome
        and c_i = the count of aligned reads for the patient at genome
        position i.
    """

    # Your code here
    patient_dict = get_fasta_dict('patients/' + patient + '.fasta')
    reference_dict = get_fasta_dict('reference_genomes/' + bacteria + '.fasta')
    ret = [0] * 15000
    mad = make_all(reference_dict.get(next(iter(reference_dict))))
    for read in patient_dict:
        r = patient_dict.get(read)
        temp = find(r, mad)
        for i in temp:
            ret[i] = ret[i] + 1
    return ret


def longest_zeros(count_vector):
    """Given a count vector, return the start and stop position of the longest
    string of zeros in the vector.

    Args:
        count_vector (list of ints): vector of aligned read counts to the
        genome. see: the return value for `read_mapper()`

    Returns:
        (tuple of (int, int)): the start and stop position of the longest
        string of zeros in the count_vector.
    """

    # Your code here

    return (-1, -1)


def align_patient_reads():
    """YOUR CODE GOES HERE..."""
    p11 = read_mapper('patient1', 'Bacteroides_ovatus')
    p2 = read_mapper()

if __name__ == '__main__':
    align_patient_reads()
