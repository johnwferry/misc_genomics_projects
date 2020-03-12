from compsci260lib import *


def solve_ultrametric_additive():
    """
    Given distance metrics, determine if they are ultrametric and/or additive.
    """

    # Distance metrics for table 1 and table 2
    dist_1 = {"1,2": 0.3, "1,3": 0.7, "1,4": 0.9,
                "2,3": 0.6, "2,4": 0.8,
                "3,4": 0.6} # fill in table 1 here

    dist_2 = {"1,2": 0.8, "1,3": 0.4, "1,4": 0.6, "1,5": 0.8,
              "2,3": 0.8, "2,4": 0.8, "2,5": 0.4,
              "3,4": 0.6, "3,5": 0.8,
              "4,5": 0.8}

    # Check if dist_1 and dist_2 are ultrametric and additive by
    # calling is_ultrametric and is_additive with the default
    # threshold value (1e-4).
    #
   # print is_ultrametric(dist_1)
    #print is_additive(dist_1)
    #print is_ultrametric(dist_2)
    #print is_additive(dist_2)
    # Your code here
    #

    # Construct the ATPA synthase distance metric table
    atpa_table = {"1,2": 0.5, "1,3": 0.5, "1,4": 0.1, "1,5": 0.4, "1,6":0.4,
                  "2,3": 0.3, "2,4": 0.5, "2,5": 0.5, "2,6": 0.5,
                  "3,4": 0.5, "3,5": 0.5, "3,6": 0.5,
                  "4,5": 0.4, "4,6": 0.4,
                  "5,6": 0.3} #  fill in ATPA synthase distance metric table

    # Determine if the ATPA synthase distance metrics
    # are ultrametric and additive
    #
    # Your code here
    #
    print (is_ultrametric(atpa_table))
    print (is_additive(atpa_table))


def is_ultrametric(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are ultrametric.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called.

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal
    Returns:
        (bool) True if the given distance metric is an ultrametric,
    False otherwise."""


    def tds(a,b):
        return dist.get(str(a) + ',' + str(b))

    """YOUR CODE GOES HERE..."""
    il = len(dist.keys())
    i = 1
    c = 0
    while c < il:
        c = c + i
        i = i + 1
    """ i is the max distance metric"""
    dl = []
    for x in range(1,i+1):
        for z in range(1,i+1):
            for y in range(1,i+1):
                if not(x == y or x == z or y == z):
                    if y > z > x:
                        xt = x
                        yt = z
                        zt = y
                    if y > x > z:
                        xt = z
                        yt = x
                        zt = y
                    if x > y > z:
                        xt = z
                        yt = y
                        zt = x
                    if x > z > y:
                        xt = y
                        yt = z
                        zt = x
                    if z > x > y:
                        xt = y
                        yt = x
                        zt = z
                    if z > y > x:
                        xt = x
                        yt = y
                        zt = z
                    if (tds(xt,zt) > max(tds(xt,yt), tds(yt,zt))) and not (is_almost_equal(tds(xt,zt), max(tds(xt,yt), tds(yt,zt)), threshold)):
                        #print ("Sequence failed at ", (xt,yt,zt), " with values: ")
                       # print ("    ",(xt,zt), tds(xt,zt))
                        #print ("    ",(xt,yt), tds(xt,yt))
                        #print ("    ",(yt,zt), tds(yt,zt))
                        return False
    return True


def is_additive(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are additive.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called.

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal

    Returns:
        (bool) Return True if the given distance metric is additive, 
        False otherwise."""

    """YOUR CODE GOES HERE..."""

    def tds(a, b):
        return dist.get(str(a) + ',' + str(b))

    """YOUR CODE GOES HERE..."""
    il = len(dist.keys())
    i = 1
    cn = 0
    while cn < il:
        cn = cn + i
        i = i + 1
    """ i is the max distance metric"""
    dl = []
    for a in range(1, i + 1):
        for b in range(a, i + 1):
            for c in range(b, i + 1):
                for d in range(c, i + 1):
                    if not(a == b or a == c or a == d or b == c or b == d or c == d):

                        if ((tds(a,b) + tds(c,d)) > max(tds(a,c) + tds(b,d), tds(a,d) + tds(b,c))) and not(is_almost_equal((tds(a,b) + tds(c,d)),max( tds(a,c) + tds(b,d), tds(a,d) + tds(b,c)), threshold)):
                            print ("Sequence failed at ", (a, b, c, d), " with values: ")
                            print ("    ", (a, b), tds(a, b))
                            print ("    ", (a, c), tds(a, c))
                            print ("    ", (a, d), tds(a, d))
                            print ("    ", (b, c), tds(b, c))
                            print ("    ", (b, d), tds(b, d))
                            print ("    ", (c, d), tds(c, d))
                            return False
    return True


def is_almost_equal(num_1, num_2, threshold=0):
    """Return true if the difference between the two parameters is negligible
    enough that the parameters can be considered equal.

    Args:
        num_1 (float/int): numeric value to compare
        num_2 (float/int): numeric value to compare
        threshold (float): maximium difference in which numeric values are 
            considered equal

    Returns:
        (bool): true if the numeric values are within the threshold
    """
    return abs(num_1 - num_2) <= threshold


if __name__ == '__main__':
    solve_ultrametric_additive()
