'''
@author: John Ferry
@date: 9/13/2019

@note: It's 6:00 am here and I started at 10:00pm last night: doing this all at once is definitely a bad idea.

'''


from orfs import *


def solve_plasmid():

    # Load the plasmid reads from fasta ensuring they are
    # in the expected order to be assembled
    #
    plasmid_dict = get_fasta_dict('plasmid.fasta')
    reads = [plasmid_dict['start']]
    plasmid_dict.pop('start')
    for l in plasmid_dict:
        reads.append(plasmid_dict[l])


    # Assemble the reads into a single-stranded linearized version of the
    # plasmid:
    plasmid_dna = simple_assembler(reads)

    # Report the length of the assembled plasmid
    #
    print (len(plasmid_dna))
    #

    # Search for ORFs in the reconstructed plasmid DNA
    #
    orfs50_plasmid = find_orfs(plasmid_dna, 50, True, True)
    print (orfs50_plasmid)
    # print summarize_orfs(orfs50_plasmid)
    ''' 
    For part h
    1. Construct a list of tuples for  (start, stop) nucleotides.
    2. Order list based on start
    3. Collapse this list by removing overlaps
    4. Calculate nucleotides by (stop-start) summation accros all tuples. 
    5. above / total = %
    
    
    (key notes for 1)
    start(does code for something so include), stop -2 (don't include itself). For +, remain in order & construct a tuple.
    For -, subtract value from the total number for pos position (and flip the order)

    '''
    t_list = []
    for orf in orfs50_plasmid: # 1
        if orf['strand'] == '+':
            first = orf['start']
            second = orf['stop']-2
            if first < 0 or second < 0 or int((second-first)/3) != (second-first)/3:
                print ("mistakes were made")
            if (first != second):
                t_list.append([first, second])
        if orf['strand'] == '-':
            first = 8541-(orf['stop']-2)
            second = 8541-(orf['start'])
            if first < 0 or second < 0 or int((second-first)/3) != (second-first)/3:
                print ("mistakes were made")
            if (first != second):
                t_list.append([first, second])

    t_list.sort(key=lambda rule: rule[0]) # 2
    #print t_list

    # mannual collapse:
    # [(285, 318), (327, 363), (375, 405), (453, 462), (579, 666), (948, 1116), (1236, 1242), (1335, 1416), (1584, 1605), (1683, 1689), (1773, 1809), (1821, 1890), (1896, 1935), (1953, 1968), (2079, 2091), (2100, 2151), (2178, 2343), (2346, 2382), (2430, 2595), (2955, 2994), (3015, 3024), (3069, 4011), (4179, 4182), (4311, 4314), (4332, 4377), (4419, 4434), (4458, 4488), (4524, 4533), (4596, 4602), (4659, 5109), (5175, 5307), (5310, 6351), (6417, 6429), (6453, 6456), (6495, 6816), (6855, 6987), (7029, 7041), (7044, 7164), (7329, 7998), (8022, 8043), (8217, 8274)]
    m_l = [(285, 318), (327, 363), (375, 405), (453, 462), (579, 666), (948, 1116), (1236, 1242), (1335, 1416), (1584, 1605), (1683, 1689), (1773, 1809), (1821, 1890), (1896, 1935), (1953, 1968), (2079, 2091), (2100, 2151), (2178, 2343), (2346, 2382), (2430, 2595), (2955, 2994), (3015, 3024), (3069, 4011), (4179, 4182), (4311, 4314), (4332, 4377), (4419, 4434), (4458, 4488), (4524, 4533), (4596, 4602), (4659, 5109), (5175, 5307), (5310, 6351), (6417, 6429), (6453, 6456), (6495, 6816), (6855, 6987), (7029, 7041), (7044, 7164), (7329, 7998), (8022, 8043), (8217, 8274)]
    f_l = []
    for tup in m_l:
        f_l.append(tup[1]-tup[0])
    # print f_l
    sum = 0
    for i in f_l:
        sum = sum + i
    #print sum
    # 5136
    # 5136/ 8541

    # for part i
    max_l = 0
    max_orf = {}
    for orf in orfs50_plasmid:
        if orf['aalength'] > max_l:
            max_l = orf['aalength']
            max_orf = orf
    dna_start = max_orf['start']
    dna_stop = max_orf['stop'] -2
    dna_seq = plasmid_dna[dna_start:dna_stop]
    prot_seq = translate(dna_seq)
    # print prot_seq
    #     PQSLRLTIRPSDPLRPLARSRSRET*RRHLFTHAT*LSDRLRWPPLTTISCRATEAGQPV*GWARGSRPILALVAQGQGTRCCQHTRCLPPRWLLRLRPILVSPQQVQLRPM*PHIIQLLRSLVRPHALALPTLELRIVLGFPQTRSRPRCPRTPI*AIITPHRPTQ*LPPLHQALPLT*LQGSSNSSLAASIHPGSRVA*LRLH*KDSIWYLRSAEASYLRKKSW*LLIRQTNHRW*RWFFCLQAADYAQKKRISRRSFDLFYGV*RSVERKLTLRDFGHEIIKKDLHLDPFKLKMKF*INLKYI*VNLV*QLPMLNQ*GTYLSDLSISFIHSCLTPRRVDNYDTGGLTIWPQCCNDTARPTLTGSR








def simple_assembler(reads):
    """Given a list of reads, as described in the problem, return an assembled
    DNA sequence, as a string. For consistency, use the first entry in the
    fragment reads list as the starting position of the returned sequence.

    For example, if we were to take in a list of three reads, 31 nucleotides
    long each. The last 15 nucleotides of each read would overlap with one
    other read, and the assembled sequence would be 48 nucleotides long with
    the sequence starting with the beginning of the first read.

    Args:
        reads (list): list of sequence strings as reads

    Returns:
         str: an assembled genomic sequence as a string starting with the first
              read in `reads`
    """

    #

    b = reads[0]
    reads.remove(b)
    while len(reads) > 0:
        for s in reads:
            if s[0:15] == b[len(b)-15:len(b)]:
                # print "adding ", len(s[15:len(s)]), " of ", len(s), " to ", len(b), " for a total ", len(b + s[15:len(s)])

                b = b + s[15:len(s)]
                reads.remove(s)
    if b[0:15] != b[len(b)-15:len(b)]:
        print ("Did an oopsie")

    return b[0:len(b)-15]


if __name__ == '__main__':
    solve_plasmid()
