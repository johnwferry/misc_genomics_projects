import sys
import random
from compsci260lib import *
from math import exp


def simulate():
    read_length = 450
    read_number = 4*(10**4)
    genome_length = 3 * (10**6)
    emp_arr = []
    nuc_arr = []
    contig_l_arr = []
    contig_c_arr = []

    iterrations = 20
    while iterrations > 0:
        print("----------------")
        print("Trial Number: ", iterrations)
        print("             ")
        iterrations = iterrations -1
        sequence_count = [0] * genome_length
        togo = read_number
        while togo > 0:
            tempindex = random.randint(0, genome_length-1)
            templength = read_length
            while templength > 0:
                sequence_count[tempindex] = sequence_count[tempindex] + 1
                tempindex = tempindex + 1
                if tempindex >= genome_length:
                    tempindex = 0
                templength = templength - 1
            togo = togo - 1
        print("Empirical Coverage: ", (sum(sequence_count)/len(sequence_count)))
        emp_arr.append((sum(sequence_count)/len(sequence_count)))
        zerocount = 0
        for i in sequence_count:
            if i == 0:
                zerocount = zerocount + 1
        print("Number of nucleotides not covered by any read: ", zerocount)
        nuc_arr.append(zerocount)

        contig_count = 0
        contig_arr = []
        contig_temp_length = 0
        cindex = 0
        while cindex < genome_length:
            if sequence_count[cindex] > 0:
                contig_temp_length = contig_temp_length + 1
            if contig_temp_length is not 0 and sequence_count[cindex] is 0:
                contig_arr.append(contig_temp_length)
                contig_count = contig_count + 1
                contig_temp_length = 0
            cindex = cindex + 1
        if contig_temp_length > 0:
            contig_arr.append(contig_temp_length)
            contig_count = contig_count +1
        print("Number of Contigs: ", contig_count)
        contig_c_arr.append(contig_count)
        print("Average Contig Length: ", (sum(contig_arr)/len(contig_arr)))
        contig_l_arr.append((sum(contig_arr)/len(contig_arr)))
    print("AVERAGE COVERAGE: ", (sum(emp_arr)/len(emp_arr)))
    print("AVERAGE NOT COVERED: ", (sum(nuc_arr)/len(nuc_arr)))
    print("AVERAGE NUM CONTIGS: ", (sum(contig_c_arr)/len(contig_c_arr)))
    print("AVERAGE AVG LEN CONTIGS: ", (sum(contig_l_arr)/len(contig_l_arr)))




if __name__ == '__main__':
    simulate()
