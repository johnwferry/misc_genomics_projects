from compsci260lib import *
from GlobalAligner import *


def solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and final value.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

    Reports:
        the optimal alignment score
        the up-most alignment achieving this score
        the down-most alignment achieving this score
    """
    dp_table = [[0] * (len(seq2)+1) for _ in range(len(seq1)+1)]
    I = len(dp_table)  # so I is 1 more than m
    J = len(dp_table[0])

    up_table = [[False] * (len(seq2)+1) for _ in range(len(seq1)+1)]
    left_table = [[False] * (len(seq2)+1) for _ in range(len(seq1)+1)]
    diag_table = [[False] * (len(seq2)+1) for _ in range(len(seq1)+1)]


    for i in range(I):
        dp_table[i][0] = -i * gap_penalty
        up_table[i][0] = True

    for j in range(J):
        dp_table[0][j] = -j * gap_penalty
        left_table[0][j] = True
    dp_table[0][0] = 0
    for i in range(1, I):
        for j in range(1, J):
            diag = subst_dict.get(seq1[i-1:i] + seq2[j-1:j]) + dp_table[i-1][j-1]
            left = dp_table[i][j - 1] - gap_penalty
            up = dp_table[i - 1][j] - gap_penalty
            if (up >= left) and (up >= diag):
                up_table[i][j] = True
                dp_table[i][j] = up
            if (diag >= left) and (diag >= up):
                diag_table[i][j] = True
                dp_table[i][j] = diag
            if (left >= diag) and (left >= up):
                left_table[i][j] = True
                dp_table[i][j] = left
    traceI = I-1
    traceJ = J-1
    topU = ""
    topB = ""
    # top traceback
    while (traceI+traceJ>0):
        if up_table[traceI][traceJ]:
            topU = seq1[traceI-1] + topU
            topB = '-' + topB
            traceI -= 1
        else:
            if diag_table[traceI][traceJ]:
                topU = seq1[traceI - 1] + topU
                topB = seq2[traceJ - 1] + topB
                traceJ -= 1
                traceI -= 1
            else:
                if left_table[traceI][traceJ]:
                    topU = '-' + topU
                    topB = seq2[traceJ - 1] + topB
                    traceJ -= 1
    traceI = I-1
    traceJ = J-1
    botU = ""
    botB = ""
    while (traceI+traceJ>0):
        if left_table[traceI][traceJ]:
            botU = '-' + botU
            botB = seq2[traceJ - 1] + botB
            traceJ -= 1
        else:
            if diag_table[traceI][traceJ]:
                botU = seq1[traceI - 1] + botU
                botB = seq2[traceJ - 1] + botB
                traceJ -= 1
                traceI -= 1
            else:
                if up_table[traceI][traceJ]:
                    botU = seq1[traceI - 1] + botU
                    botB = '-' + botB
                    traceI -= 1
    print("Max Score: ", dp_table[I-1][J-1])
    print(" ")
    print(topU)
    print(topB)
    print(" ")
    print(botU)
    print(botB)










if __name__ == '__main__':
    tempdict = get_fasta_dict('atpa_Hs.fasta')

    seq1 = tempdict.values()[0]
    tempdict = get_fasta_dict('atpa_Ec.fasta')

    seq2 = tempdict.values()[0]
    match = 2
    mismatch = -1
    gap_penalty = 8
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

    solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty)
