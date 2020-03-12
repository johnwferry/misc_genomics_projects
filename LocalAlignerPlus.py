from compsci260lib import *
from GlobalAligner import *


def solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
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
    zero_table = [[False] * (len(seq2)+1) for _ in range(len(seq1)+1)]

    for i in range(I):
        up = -i * gap_penalty
        zero = 0
        if (zero >= up):
            dp_table[i][0] = zero
            zero_table[i][0] = True
        if (up >= zero):
            dp_table[i][0] = up
            up_table[i][0] = True

    for j in range(J):
        left = -j * gap_penalty
        zero = 0
        if (zero >= up):
            dp_table[0][j] = zero
            zero_table[0][j] = True
        if (left >= zero):
            dp_table[0][j] = left
            zero_table[0][j] = True
    dp_table[0][0] = 0
    max_value = 0
    max_cell = [[0, 0]]
    for i in range(1, I):
        for j in range(1, J):
            diag = subst_dict.get(seq1[i-1:i] + seq2[j-1:j]) + max(dp_table[i-1][j-1], 0)
            left = dp_table[i][j - 1] - gap_penalty
            up = dp_table[i - 1][j] - gap_penalty
            if (up >= left) and (up >= diag and (up >= 0)):
                up_table[i][j] = True
                dp_table[i][j] = up
            if (diag >= left) and (diag >= up) and (diag >= 0):
                diag_table[i][j] = True
                dp_table[i][j] = diag
            if (left >= diag) and (left >= up) and (left >=0):
                left_table[i][j] = True
                dp_table[i][j] = left
            if (0 >= left) and (0 >= diag) and (0 >= up):
                zero_table[i][j] = True
                dp_table[i][j] = 0

            if dp_table[i][j] == max_value:
                max_cell.append([i,j])
            if dp_table[i][j] > max_value:
                max_cell = [[i,j]]
                max_value = dp_table[i][j]
    print("Max Score: ", max_value)
    print("The following are all max tracebacks ", max_cell)
    traceI = max_cell[0][0]
    traceJ = max_cell[0][1]
    topU = ""
    topB = ""
    topF = max_cell[0]
    topF = max_cell[0]
    topS = [0,0]
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
                else:
                    if zero_table[traceI][traceJ]:
                        topS = [traceI,traceJ]
                        traceJ = 0
                        traceI = 0

    bot_cell = max_cell[len(max_cell)-1]
    traceI = bot_cell[0]
    traceJ = bot_cell[1]
    botU = ""
    botB = ""
    botF = bot_cell
    botS = [0, 0]
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
                else:
                    if zero_table[traceI][traceJ]:
                        botS = [traceI,traceJ]
                        traceJ = 0
                        traceI = 0


    print(" ")
    print("   Starts at ", topS)
    print(topU)
    print(topB)
    print("   Ends at ", topF)
    print(" ")
    print("   Starts at ", botS)
    print(botU)
    print(botB)
    print("   Ends at ", botF)



if __name__ == '__main__':
    seq1 = "MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATPEVVSKIAQYKRECPSIFAWEIRDRLLSEGVCTNDNIPSVSSINRVLRNLASEKQQMGADGMYDKLRMLNGQTGSWGTRPGWYPGTSVPGQPTQDGCQQQEGGGENTNSISSNGEDSDEAQMRLQLKRKLQRNRTSFTQEQIEALEKEFERTHYPDVFARERLAAKIDLPEARIQVWFSNRRAKWRREEKLRNQRRQASNTPSHIPISSSFSTSVYQPIPQPTTPVSSFTSGSMLGRTDTALTNTYSALPPMPSFTMANNLPMQPPVPSQTSSYSCMLPTSPSVNGRSYDTYTPPHMQTHMNSQPMGTSGTTSTGLISPGVSVPVQVPGSEPDMSQYWPRLQ"
    seq2 = "MRNLPCLGTAGGSGLGGIAGKPSPTMEAVEASTASHRHSTSSYFATTYYHLTDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATAEVVSKISQYKRECPSIFAWEIRDRLLQENVCTNDNIPSVSSINRVLRNLAAQKEQQSTGSGSSSTSAGNSISAKVSVSIGGNVSNVASGSRGTLSSSTDLMQTATPLNSSESGGASNSGEGSEQEAIYEKLRLLNTQHAAGPGPLEPARAAPLVGQSPNHLGTRSSHPQLVHGNHQALQQHQQQSWPPRHYSGSWYPTSLSEIPISSAPNIASVTAYASGPSLAHSLSPPNDIESLASIGHQRNCPVATEDIHLKKELDGHQSDETGSGEGENSNGGASNIGNTEDDQARLILKRKLQRNRTSFTNDQIDSLEKEFERTHYPDVFARERLAGKIGLPEARIQVWFSNRRAKWRREEKLRNQRRTPNSTGASATSSSTSATASLTDSPNSLSACSSLLSGSAGGPSVSTINGLSSPSTLSTNVNAPTLGAGIDSSESPTPIPHIRPSCTSDNDNGRQSEDCRRVCSPCPLGVGGHQNTHHIQSNGHAQGHALVPAISPRLNFNSGSFGAMYSNMHHTALSMSDSYGAVTPIPSFNHSAVGPLAPPSPIPQQGDLTPSSLYPCHMTLRPPPMAPAHHHIVPGDGGRPAGVGLGSGQSANLGASCSGSGYEVLSAYALPPPPMASSSAADSSFSAASSASANVTPHHTIAQESCPSPCSSASHFGVAHSSGFSSDPISPAVSSYAHMSYNYASSANTMTPSSASGTSAHVAPGKQQFFASCFYSPWV"
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

    solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty)
