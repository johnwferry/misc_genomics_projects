'''
@author: John Ferry jwf17
@date: 9/12/2019

@note:

'''


from compsci260lib import *   # this will also import the sys and re modules


def solve_cloning():
    """Function demonstrating the procedure of extracting a gene and inserting
    into a plasmid using restriction enzymes."""

    aim2_fasta_path = 'aim2_plus_minus_1kb.fasta'
    pRS304_fasta_path = 'pRS304.fasta'

    # Read in the aim2 genomic sequence from the fasta file along with its
    # upstream and downstream regions.
    aim2_genomic = get_fasta_dict(aim2_fasta_path)["aim2"]
    # Store the beginning and the end of the Aim2 gene as Python indices.
    aim2_beg = 1001
    aim2_end = 1741
    print("Begin:", aim2_beg)
    print("End: ", aim2_end)
    # calc nuc length & amino acids
    aim2_length_nuc = aim2_end - aim2_beg
    print("Nuc. Length: ", aim2_length_nuc)
    aim2_length_amino = aim2_length_nuc / 3
    print("Amino Length: ", aim2_length_amino)

    # Define regular expression terms for each restriction enzyme
    r_enzymes = get_restriction_enzymes_regex()

    #  print r_enzymes

    # Store coordinates of restriction sites found upstream, downstream, and
    # within the aim2 gene
    r_enzyme_sites = find_aim2_restriction_enzymes(
        aim2_beg, aim2_end, aim2_genomic)

    # Report the found restriction enzyme sites
    #
    print("Aim2 Enzyme sites: ", r_enzyme_sites)
    #

    # Input the pRS304 plasmid selecting the MCS
    # Store the MCS in a new variable.
    #
    temp_mcs = get_fasta_dict(pRS304_fasta_path)

    mcs_start = 1887
    mcs_end = 1994
    prs304_mcs = temp_mcs["pRS304"][mcs_start:mcs_end]

    # Find and report the restriction enzymes in the MCS
    p_enzyme_sites = find_pRS304_MCS_restriction_sites(prs304_mcs, mcs_start)
    #
    print("MCS Enzyme sites: ", p_enzyme_sites)
    #

    # Extract aim2 gene and insert into the plasmid, report the length
    #
    #
    aim2_extract = aim2_genomic[688:2361]
    print("Beginning extract nuc pos ", 688)
    print("Ending Extract nuc pos", 2361)
    #   print "aim ", len(aim2_extract)
    plasmid_int = temp_mcs["pRS304"][0:1967] + (aim2_extract)
    #  print "raw plas", len(temp_mcs["pRS304"])
    plasmid_int = plasmid_int + temp_mcs["pRS304"][1972:4267]
    print("Plasmid w/ aim2 length ", len(plasmid_int))
    #


def get_restriction_enzymes_regex():
    """Returns a dictionary of restriction enzyme regular expressions for
    searching in genomic sequences.

    This function should be used for find_aim2_restriction_enzymes and
    find_pRS304_MCS_restriction_sites.
    """

    r_enzymes = {
        "BamHI": ["GGATCC"],
        "BstYI": ["AGATCC", "AGATCT", "GGATCC", "GGATCT"],
        "SalI": ["GTCGAC"],
        "SpeI": ["ACTAGT"],
        "SphI": ["GCATGC"],
        "StyI": ["CCAAGG", "CCATGG", "CCTAGG", "CCTTGG"],
    }
    return r_enzymes


def find_aim2_restriction_enzymes(aim2_beg, aim2_end, aim2_genomic):
    """Finds the restriction enzyme sites in the aim2 genomic sequence. Stored
    as a dictionary of lists. Each restriction enzyme key corresponds to a list
    of dictionaries containing relevant information for a found site.

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched in the genome
        start (int): the start nucleotide position in the genome
        end (int): the ending position
        location (str): the position of the site relative the aim2 gene.
            (must be "upstream", "downstream", or "within")

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC',
                'location': 'upstream'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC',
                'location': 'downstream'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC',
                'location': 'within'
            }
        ]
    }
    """

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    r_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}

    # var for later use
    aim2_s = aim2_genomic
    e = len(aim2_s)-7
    if e < 0:
        e = 0
    # hardcode 6 as recognition site length
    for i in range(0, e):  # check every index in the entire aim2 length
        for z in r_enzymes:  # check every enzyme in the dict
            for c in r_enzymes[z]:  # check every cut sequence in enzyme
                #  print "Checking ", c.lower(), " and ", aim2_s[i:i+6]
                if c.lower() == aim2_s[i:i+6]: # if cut matches string
                    #  print "Passed ", c.lower(), " and ", aim2_s[i:i+6]
                    if i+6 < aim2_beg:  # upstream
                        r_enzyme_sites[z].append({'start': i, 'end': i+5, 'sequence': c, 'location': "upstream"})
                    if i > aim2_end:  # downstream
                        r_enzyme_sites[z].append({'start': i, 'end': i+5, 'sequence': c, 'location': "downstream"})
                    if i <= aim2_end and i+6 >= aim2_beg:  # within
                        r_enzyme_sites[z].append({'start': i, 'end': i+5, 'sequence': c, 'location': "within"})
    return r_enzyme_sites


def find_pRS304_MCS_restriction_sites(prs304_mcs, mcs_start):
    """Finds the restriction sites for for the MCS of pRS304. Stored as a
    dictionary of lists. Each restriction enzyme key corresponds to a list of
    dictionaries containing relevant information for a found site in the MCS.

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched in the MCS
        start (int): the start nucleotide position in the MCS
        end (int): the ending position in the MCS

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC'
            }
        ]
    }
    """

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    p_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}
    # var for later use
    e = len(prs304_mcs) - 7
    # hardcode 6 as recognition site length
    for i in range(0, e):  # check every index in the entire length
        for z in r_enzymes:  # check every enzyme in the dict
            for c in r_enzymes[z]:  # check every cut sequence in enzyme
                #  print "Checking ", c.lower(), " and ", prs304_mcs[i:i+6]
                if c.lower() == prs304_mcs[i:i + 5]:  # if cut matches string
                    #  print "Passed ", c.lower(), " and ", prs304_mcs[i:i+6]
                    p_enzyme_sites[z].append({'start': i + mcs_start, 'end': i + 6 + mcs_start, 'sequence': c})

    return p_enzyme_sites


if __name__ == '__main__':
    solve_cloning()
