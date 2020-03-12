'''
@author: John Ferry
@date: jwf17

@note:

'''


from compsci260lib import *   # this will also import the sys and re modules


def solve_orfs():
    """Your code goes here..."""
    data = get_fasta_dict('sars.fasta')["sars"]
    orfs10 = find_orfs(data, 10)
    orfs50 = find_orfs(data, 50)
    orfs70 = find_orfs(data, 70)
    s_orfs10 = summarize_orfs(orfs10)
    s_orfs50 = summarize_orfs(orfs50)
    s_orfs70 = summarize_orfs(orfs70)
    print ("min10aa Number: ", s_orfs10[0])
    print ("min10aa Avg: ", s_orfs10[1])
    print ("min50aa Number: ", s_orfs50[0])
    print ("min50aa Avg: ", s_orfs50[1])
    print ("min70aa Number: ", s_orfs70[0])
    print ("min70aa Avg: ", s_orfs70[1])


def summarize_orfs(orfs):
    """Summarize ORFs identified from the find_orfs procedure as a count of the
    number of found orfs and the average length of the found ORFs (in amino
    acids)

    Args:
        orfs (list): a list of dictionaries of found ORFs

    Returns:
         tuple: (The number of ORFs found (int), Average ORF length (float))
    """

    #
    lengths = []

    for orf in orfs: #
        lengths.append(orf['aalength'])
    sum = 0
    avg = 0
    for x in lengths: # avging
        sum += x
        avg = sum / len(lengths)

    return len(lengths), avg  # replaced with the tuple described above


def find_orfs(seq, min_length_aa=0, is_double_stranded=False,
              is_circular=False):
    """This is a function for finding sufficiently long ORFs in all reading
    frames in a sequence of DNA or RNA.  By default, the sequence is assumed
    to be single-stranded, but if you call this function with
    is_double_stranded = True, it will look for ORFs on both strands.

    The function takes as input parameters: a string representing a genomic
    sequence, the minimum length (in amino acids) for an ORF before it will be
    returned (which defaults to 0), whether the sequence should be considered
    double-stranded, and whether the sequence is circular (like a plasmid).

    The last two parameters default to False; so if default values are used,
    the function will assume a single-stranded linear sequence.

    Args:
        seq (str): a genomic sequence
        min_length_aa (int): minimum length of found ORFs in amino acids
        is_double_stranded (bool): whether the sequence should be considered
                                   double-stranded
        is_circular (bool): whether the sequence is circular (like a plasmid)

    Returns:
        list: of dictionaries with information on each ORF found.

    Where each ORF found is represented by a dictionary with
    the following keys:
        frame (int): The nucleotide offset in which the ORF was found. (Must be
        0, 1, or 2)
        stop (int): the nucleotide position of the end of the ORF
        start (int): the nucleotide position of the start of the ORF
        stopcodon (str): the nucleotide triplet of the stop codon
        nlength (int): the length (in nucleotides) of the ORF
        strand (str): the strand of the found ORF (Must be '+' or '-')

    A valid return list may look something like this:
    [
        {
            'frame': 0,
            'stop': 13413,
            'aalength': 4382,
            'start': 265,
            'stopcodon': 'UAA',
            'nlength': 13149,
            'strand': '+'
        },
        {
            'frame': 0,
            'stop': 27063,
            'aalength': 221,
            'start': 26398,
            'stopcodon': 'UAA',
            'nlength': 666,
            'strand': '-'
        }
    ]
    """

    #
    ret = [] # return var
    if not is_circular:
        s = seq.lower()
        f = 0 # frame var
        start_c = 'aug'
        end_c = {'uag', 'uga', 'uaa'}
        i = 0
        while i < len(s)-2: # while codons remain in s
            if s[i:i+3] == start_c: # if this codon is start, move deeper
                x = i # set secondary index to current
                while x < len(s)-2: # while codons remain in s
                    for c in end_c: # check all stop codons
                        if s[x:x+3] == c and (x-i-3)/3 >= min_length_aa: # if it is a stop codon and aa length is satisfy
                            ret.append({'frame': f, 'stop': x+2, 'aalength': int((x-i-3)/3), 'start': i, 'stopcodon': c, 'nlength': (x-i-3), 'strand': '+'})
                            i = x+2 # i gets iterrated later, move it up to stop codon so it moves forward
                            x = len(s)
                            f -= 1 # f gets iterrated later, this corrects the frame mis-allignment
                    x += 3 # iterrate by 3 so that you can remain in frame set
            i += 1 # move along 1 nuc
            if f < 2: # iterrates f looping it back to 0
                f += 1
            else:
                f = 0
        return ret
    if is_circular and is_double_stranded:
        l = int(len(seq)/3)
        start_c = 'atg'
        end_c = {'tag', 'tga', 'taa'}
        for f in range(0,2):
            if f == 0:
                s = seq
            else:
                s = seq[f:len(seq)] + seq[0:f]
            i = 0
            passed_initial_start_c = False
            initial_start_c_num = -1
            while not passed_initial_start_c:
                c = s[0:3]
                s = s[3:len(s)] + s[0:3]
                if int(len(s)/3) != l:
                    print ("ya dun did an oopsie")
                i = i + 1
                if c == start_c:
                    f_end = False
                    if initial_start_c_num == -1:
                        initial_start_c_num = i
                    if i >= l + initial_start_c_num or i >= 3*l:
                        passed_initial_start_c = True
                    if passed_initial_start_c:
                        break
                    x = i
                    while not f_end:
                        c = s[0:3]
                        s = s[3:len(s)] + s[0:3]
                        i = i + 1
                        for end in end_c:
                            if end == c:
                                ret.append({
                                    'frame': f,
                                    'stop': 3*i -1,  # stop at last part of stop codon
                                    'aalength': i - x -1,  # length is ignoring start and stop codons
                                    'start': 3*(x -1),
                                    'stopcodon': end,
                                    'nlength': 3 *(i - x - 1),  # length is ignoring start and stop codons
                                    'strand': '+'
                                })
                                f_end = True
        r_seq = seq[::-1] # reversed seq for flip side
        n_seq = '' # flip all the atcg's to their compliments to complete the reverse 'negated' strand
        for cha in r_seq:
            if cha == 'a':
                n_seq += 't'
            if cha == 't':
                n_seq += 'a'
            if cha == 'g':
                n_seq += 'c'
            if cha == 'c':
                n_seq += 'g'
        if len(n_seq) != len(seq):
            print (" reverse oopsie ")
        for f in range(0,2): # run like normal
            if f == 0:
                s = n_seq
            else:
                s = n_seq[f:len(seq)] + n_seq[0:f]
            i = 0
            passed_initial_start_c = False
            initial_start_c_num = -1
            while not passed_initial_start_c:
                c = s[0:3]
                s = s[3:len(s)] + s[0:3]
                if int(len(s)/3) != l:
                    print ("ya dun did an oopsie")
                i = i + 1
                if c == start_c:
                    f_end = False
                    if initial_start_c_num == -1:
                        initial_start_c_num = i
                    if i >= l + initial_start_c_num or i >= 3*l:
                        passed_initial_start_c = True
                    if passed_initial_start_c:
                        break
                    x = i
                    while not f_end:
                        c = s[0:3]
                        s = s[3:len(s)] + s[0:3]
                        i = i + 1
                        for end in end_c:
                            if end == c:
                                ret.append({
                                    'frame': f,
                                    'stop': 3*i -1,  # stop at last part of stop codon
                                    'aalength': i - x -1,  # length is ignoring start and stop codons
                                    'start': 3*(x -1),
                                    'stopcodon': end,
                                    'nlength': 3 *(i - x - 1),  # length is ignoring start and stop codons
                                    'strand': '-'
                                })
                                f_end = True
    return ret








if __name__ == '__main__':
    solve_orfs()
