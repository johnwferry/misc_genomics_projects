'''
@author:
@date:

@note:

'''

from compsci260lib import *





def solve_pebble(grid_file):
    """Code for the "Pebble Beach" problem. This problem involves implementing
    an O(n) dynamic programming algorithm for computing the maximum value of
    the placement of pebbles under the constraint that no pebbles can be
    vertically or horizontally adjacent.

    Args: grid_file (str): a string with the name of the file that contains
          the grid of values. Each line of that file should contain a row
          of four integers, separated by tabs

    Returns: the maximal score for the optimal pebble placements
    """
    configs = [
        [1, 1, 0, 0, 1],
        [2, 1, 0, 1, 0],
        [3, 0, 1, 0, 1],
        [4, 1, 0, 0, 0],
        [5, 0, 1, 0, 0],
        [6, 0, 0, 1, 0],
        [7, 0, 0, 0, 1],
        [8, 0, 0, 0, 0]
    ]
    # is b compatible with a

    def compatible(a, b):
        if a[0] is 1:
            return not((b[0] is 1) or (b[0] is 2) or (b[0] is 3) or (b[0] is 4) or (b[0] is 7))
        if a[0] is 2:
            return not((b[0] is 2) or (b[0] is 1) or (b[0] is 4) or (b[0] is 6))
        if a[0] is 3:
            return not((b[0] is 3) or (b[0] is 1) or (b[0] is 5) or (b[0] is 7))
        if a[0] is 4:
            return not((b[0] is 4) or (b[0] is 1) or (b[0] is 2))
        if a[0] is 5:
            return not((b[0] is 5) or (b[0] is 3))
        if a[0] is 6:
            return not((b[0] is 6) or (b[0] is 2))
        if a[0] is 7:
            return not((b[0] is 7) or (b[0] is 3) or (b[0] is 1))
        if a[0] is 8:
            return True
        #  print "shit isn't supposed to reach here"
        return False



    def compute_value(rw):

        def multrw(cr):
            mwret = [0, 0, 0, 0, 0]
            mwret[0] = cr[0]
            for mwi in range(0, 4):
                mwret[mwi+1] = cr[mwi+1] * rw[mwi]
            return mwret
        cvret = []
        for cvr in configs:
            cvret.append(multrw(cvr))

        return cvret



    """Your code goes here..."""
    files = open(grid_file, "r")
    line = "start"
    prev = [
        [1, 0, 0, 0, 0],
        [2, 0, 0, 0, 0],
        [3, 0, 0, 0, 0],

        [4, 0, 0, 0, 0],
        [5, 0, 0, 0, 0],
        [6, 0, 0, 0, 0],
        [7, 0, 0, 0, 0],

        [8, 0, 0, 0, 0]
    ]
    scorel = [0, 0, 0, 0, 0, 0, 0, 0]
    while line is not "":
        line = files.readline()
        row = list(map(int, line.split()))
        if len(row) is 0:
            break

        cur = compute_value(row)

        def calc_score(a, b, pr):
            ascore = sum(a[1:len(a)])
            if not compatible(a, pr):
                return max(ascore, b)
            return ascore + b
        templist = [0, 0, 0, 0, 0, 0, 0, 0]
        tempscore = [0, 0, 0, 0, 0, 0, 0, 0]
        for icalc in range(0, 8):
            for jcalc in range(0, 8):
                templist[jcalc] = calc_score(cur[icalc], scorel[jcalc], prev[jcalc])

            tempscore[icalc] = max(templist)
            templist = [0, 0, 0, 0, 0, 0, 0, 0]
        for repl in range(0, 8):
            scorel[repl] = tempscore[repl]
        prev = cur












    #  Return the maximum value of the placement of pebbles
    return max(scorel)


if __name__ == '__main__':
    max_score = solve_pebble('grid.txt')
    print("The max score for this grid is %d" % max_score)
