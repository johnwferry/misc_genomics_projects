import sys
from compsci260lib import *


def solve_global_aligner(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and final value.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): score of inserting a gap
    """

    # create the DP table
    dp_table = [[0] * (len(seq2)+1) for _ in range(len(seq1)+1)]

    # Compute the score of the optimal global alignment
    max_value = global_aligner(subst_dict, gap_penalty, dp_table, seq1, seq2)

    # Display the dp table
    display_dp_table(seq1, seq2, dp_table)

    print("\n\nThe score of an optimal global alignment is " + str(max_value))


def create_subst_matrix_dna(match, mismatch):
    """Create a substitution matrix for DNA sequences. This method
    should be used with create_subst_matrix_dict().

    Args:
        match (int): score for a match
        mismatch (int): score for a mismatch

    Returns:
        (list of strings): representing the substitution matrix for DNA. Scores
        are separated by a single space.
    """
    subst_matrix = []
    subst_matrix.append("A C G T")
    subst_matrix.append("%s %s %s %s" % (match, mismatch, mismatch, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, match, mismatch, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, mismatch, match, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, mismatch, mismatch, match))

    return subst_matrix


def create_subst_matrix_aa(blosum_path):
    """Create a substitution matrix for amino acid sequences. This method
    should be used with create_subst_matrix_dict().

        Args:
            blosum_path (str): path to the BLOSUM file to load

        Returns:
            (list of strings): representing the substitution matrix for amino
            acids. Scores are separated by spaces.
        """

    subst_matrix = load_matrix(blosum_path)
    return subst_matrix


def global_aligner(subst_dict, gap_penalty, dp_table, seq1, seq2):
    """A dynamic programming algorithm that takes two sequences and returns the
    score of the optimal alignment.

    Args:
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """
    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)      # so I is 1 more than m
    J = len(dp_table[0])   # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    for i in range(I):
        dp_table[i][0] = -i * gap_penalty

    for j in range(J):
        dp_table[0][j] = -j * gap_penalty
    dp_table[0][0] = 0
    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    for i in range(1, I):
        for j in range(1, J):

            dp_table[i][j] = max(subst_dict.get(seq1[i-1:i] + seq2[j-1:j]) + dp_table[i-1][j-1],
                                  dp_table[i][j - 1] - gap_penalty,
                                  dp_table[i - 1][j] - gap_penalty
            )



    """YOUR CODE GOES HERE..."""

    # The optimal score is found at the lower right corner of the dp table:
    return dp_table[I-1][J-1]


def compute_global_aligner_score_linspace(subst_dict, gap_penalty, seq1, seq2):
    """Extra Challenge: compute the score (not the actual alignment) of the
       best global alignment in O(mn) time and only O(m) or O(n) space

    Args:
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score

    """

    """YOUR CODE GOES HERE..."""

    return -1


########################################################################
########################################################################
# Below are some functions that may be handy
########################################################################
########################################################################


def display_dp_table(seq1, seq2, dp_table):
    """A function that displays the dp table once it's been computed."""
    # Add gap character in front of the two sequences
    seq1 = "-" + seq1
    seq2 = "-" + seq2

    # Trying to get the columns lined up; this works as long as |values| < 99,
    # and can be modified for bigger values
    print("\n\n\t\t ")

    for i in range(len(seq2)):
        print(seq2[i] + "\t " + ' ')

    for i in range(len(seq1)):
        print("\n\t" + seq1[i] + ' ')

        for j in range(len(seq2)):
            adjustment = '  '
            if dp_table[i][j] > 9:
                adjustment = ' '
            if dp_table[i][j] < 0:
                adjustment = ' '
            if dp_table[i][j] < -9:
                adjustment = ''
            print("\t", adjustment, dp_table[i][j], ' ')


def is_dna(in_str):
    """Determine if string is a valid DNA sequence; it's valid if it only
    contains {ACGT} characters."""
    if re.findall(r"^[ACGT]+$", in_str.upper()):
        return True
    else:
        return False


def is_protein(in_str):
    """Determine if string is a valid protein sequence; it's valid if it only
    contains {ACGTRNDQEHILKMFPSWYV} characters."""
    if re.findall(r"^[ACGTRNDQEHILKMFPSWYV]+$", in_str):
        return True
    else:
        return False


def validate_sequences(seq1, seq2):
    """Determines if a given pair of sequences are both valid DNA sequences,
    both valid protein sequences, or an invalid pair of sequences."""
    if is_dna(seq1) and is_dna(seq2):
        return 1
    elif is_protein(seq1) and is_protein(seq2):
        return 2
    else:
        return 0


def get_sequence():
    """A console interface for sequence input.

    Returns a string sequence the user entered or from the file the user
    specified.
    """
    print("\nChoose one of the following options")
    print("\n\t1 The sequence will be entered by hand")
    print("\n\t2 The sequence will be read from a file")
    print("\nEnter choice: ")

    while True:
        try:
            choice = int(input())
            break
        except Exception:
            print("Please give valid choice: ")

    seq = ""
    if choice == 1:
        print("\nEnter a sequence: ")
        sys.stdout.flush()
        seq = input()
        seq = seq.strip().upper()
    else:
        print("\nEnter the filename: ")
        sys.stdout.flush()
        file_name = input()
        file_name = file_name.strip()
        seq = list(get_fasta_dict(file_name).values())[0]

    return seq


def create_subst_matrix_dict(rows):
    """Reads in a list of rows of strings.

    The first row labels the character types (nucleotide types if DNA,
    amino acid types if protein). The remaining rows are strings that
    contain the values in the substitution matrix.  Function transforms
    these rows of strings into a matrix, and then transforms that into a
    dictionary where the key is a pair of characters (literally, a
    string of length two) and the value is the score of aligning those
    characters. The reason for doing this is so that we can index into a
    data structure using the characters themselves.  Function then
    returns the dictionary.  Gap score is not included in any of this;
    that's handled separately.

    Args:
        (list of strings): representation of the substitution matrix

    Returns:
        (dictionary string -> int): dictionary representation of the
        substitution matrix mapping a string representing the row and column
        characters to an alignment score of the two characters:

            <first character> + <second character> -> <alignment score>

        The returned dictionary may look something like this:
        {
            'AA': 1,
            'AT': 0,
            ...
        }

    """
    matrix = []

    # Remove whitespace in front of the first row
    rows[0].strip()

    # The first row of the matrix contains a sequence of characters.  Store
    # them as a list. These characters are labels for the rows & columns of the
    # substitution matrix.
    reference_list = rows[0].split()

    # Now get rid of the first row
    rows.pop(0)

    for row in rows:
        # Again, remove whitespace in front of rows
        row = row.strip()
        cols = row.split()
        matrix.append(cols)

    # Create the dictionary
    subst_dict = {}
    for i in range(len(reference_list)):
        for j in range(len(reference_list)):
            # In the following line, the + is used to join the two characters
            # into a string of length 2.
            subst_dict[reference_list[i] +
                       reference_list[j]] = int(matrix[i][j])

    return subst_dict


def load_matrix(in_file):
    """Function to read a substitution matrix in from a file."""

    with open(in_file, 'r') as f:

        # Check to see if the file exists
        while f is None:
            print("File does not exist. \nEnter the file name again: ")
            in_file = input()
            in_file = in_file.strip()
            f = open(in_file, 'r')

        lines = f.readlines()

    return lines


if __name__ == '__main__':
    # Please update the input parameters according to the problem description.
    # The sequences can either be sequence strings or you can read them from
    # the fasta file.
    seq1 = "GAATCGGA"
    seq2 = "TAGTA"
    match = 2
    mismatch = -1
    gap_penalty = 2
    seq_type = validate_sequences(seq1, seq2)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit('Input sequences are of different types')

    # Obtain a dictionary of scores for aligning a pair of characters
    subst_dict = create_subst_matrix_dict(subst_matrix)

    solve_global_aligner(seq1, seq2, subst_dict, gap_penalty)
