from compsci260lib import *

# Note that the '$' character will be used to designate the end of a
# given string.


def forward_bwt(seq):
    """forward_bwt(seq) takes as input a string containing the EOF character to
    which the BWT must be applied. The method should then return the result of
    the BWT on the input string.

    For example:
        forward_bwt('GATTACA$') --> 'ACTGA$TA'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """
    bwm = []
    endi = len(seq)
    tempseq = seq
    i = 0
    while i < endi:
        tempseq = tempseq[1: endi] + tempseq[0:1]
        bwm.append(tempseq)
        i = i + 1
    bwm = sorted(bwm)
    ret = ""
    for s in bwm:
        ret = ret + s[endi-1:endi]
    return ret


def reverse_bwt(seq):
    """reverse_bwt(seq) takes as input a string containing the EOF character to
    which the reverse of the BWT must be applied. The method should then return
    the result of the reversal on the input string.

    For example:
        reverse_bwt('ACTGA$TA') --> 'GATTACA$'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """
    slen = len(seq)
    def reverse_bwt_inner(bwm):
        count = 2
        while count < slen:
            if len(bwm[0]) is slen:
                return bwm
            for i in range(0, slen):
                bwm[i] = seq[i] + bwm[i]
            bwm = sorted(bwm)
            count = count + 1
        return bwm
    tempbmw = []
    for s in sorted(seq):
        tempbmw.append(s)
    for ind in range(0, slen):
        tempbmw[ind] = seq[ind] + tempbmw[ind]
    tempbmw = sorted(tempbmw)
    solvedbwm = reverse_bwt_inner(tempbmw)

    return solvedbwm[0][1:slen] + '$'


def solve_bwt():

    # example sequences for forward_bwt() and reverse_bwt()
    # you may change them to different sequences to test your code.

    seq1 = 'GATTACA$'
    seq2 = 'ACTGA$TA'
    forward_bwt(seq1)
    reverse_bwt(seq2)
    seq3 = 'CGGACTAACGGACTAACGGACTAACGGACTAA$'
    print(forward_bwt(seq3))

    fp = open('mystery.txt', 'r', encoding='utf-8')
    seq4 = fp.read()
    print(reverse_bwt(seq4))
    fp.close()



if __name__ == '__main__':
    solve_bwt()



