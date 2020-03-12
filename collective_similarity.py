'''
@author: jwf17
@date: 9/27/19

@note:

'''

from compsci260lib import *
import timeit
import random


def brute_force(score_list):
    """Get the maximum similarity score using brute force.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """

    max_s = 0
    for i in range(len(score_list)+1):
        for j in range(i+1, len(score_list) + 1):
            temp = score_list[i: j]
            if sum(temp) > max_s:
                max_s = sum(temp)
    return max_s



def divide_conquer(score_list):
    """Get the maximum similarity score using divide and conquer.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """

    # compress the list to make it an easier problem
    compressed_list = []
    previous = 0
    for i in range(len(score_list)):
        snapshot = previous
        if snapshot >= 0:
            if score_list[i] >= 0:
                previous = previous + score_list[i]
            if score_list[i] < 0:
                compressed_list.append(previous)
                previous = score_list[i]
        if snapshot < 0:
            if score_list[i] >= 0:
                compressed_list.append(previous)
                previous = score_list[i]
            if score_list[i] < 0:
                previous = previous + score_list[i]



    # returns the max sum given the list and whether that sum is complete in a tuple
    # [sum, completeLeft, CompleteRIght, use both]
    # completeness is defined as being able to combine with the left or right node

    def divide_sum(s_l):

        if len(s_l) == 0:
            return [0, True, True, 0]
        if len(s_l) == 1:
            return [s_l[0], True, True, 0]

        def find_first_negative(sl):
            for j in range(len(sl)):
                if sl[j] < 0:
                    return j
            return -1

        index = find_first_negative(s_l)
        if index == -1:
            return [sum(s_l), True, True, 0]

        # recursion
        negative = s_l[index]
        left = s_l[0: index]
        right = s_l[index+1: len(s_l)]
        div_left = divide_sum(left)
        div_right = divide_sum(right)

        if div_left[3] < 0 or div_right[3] < 0:
            return [brute_force(s_l), False, False, 0]

        if div_left[2] and div_right[1]:
            default_max = div_right[0] + div_left[0] + negative
            if max([div_right[0], div_left[0], default_max]) == default_max:
                return [default_max, div_left[1], div_right[2], 0]
        if div_right[0] > div_left[0]:
            return [div_right[0], False, div_right[2], 0]
        if div_left[0] > div_right[0]:
            return [div_left[0],  div_left[1], False, 0]
        return [div_right[0], div_left[1], div_right[2], negative]  # FLAG for equal
    result = divide_sum(compressed_list)
    return max([result[0], 0])





def linear(score_list):
    """Get the maximum similarity score in linear time.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """
    c = score_list
    MAX_INCLUDING_HERE = 0
    MAX_SO_FAR = 0
    for i in range(len(c)):
        MAX_INCLUDING_HERE = MAX_INCLUDING_HERE + c[i]
        if MAX_INCLUDING_HERE < 0:
            MAX_INCLUDING_HERE = 0
        else:
            if MAX_SO_FAR <= MAX_INCLUDING_HERE:
                MAX_SO_FAR = MAX_INCLUDING_HERE

    return MAX_SO_FAR  # replace with the computed maximal score


if __name__ == '__main__':
    """You can use this to test the correctness of your code by using
    sample_list as an input to each function."""

    sample_list = [2, -3, -4, 4, 8, -2, -1, 1, 10, -5]
    brute_force(sample_list)
    divide_conquer(sample_list)
    linear(sample_list)

    """ This part below is used to test the runtime of your code, an example is
    given below for brute force algorithm with a random list of length 100.
    You will have to measure the runtime of each algorithm on every input size
    given in the problem set. """

    """
    allowed_scores = [i for i in range(-10,11)]
    random_list = [random.choice(allowed_scores) for x in range(100)]
    bruteforce_runtime = timeit.timeit('brute_force(random_list)',
        setup="from __main__ import brute_force, random_list", number=1)
    """
    """
    temp = []
    y = 5
    allowed_scores = [i for i in range(-10, 11)]
    random_list = [random.choice(allowed_scores) for x in range(10**y)]
    bruteforce_runtime = timeit.timeit('brute_force(random_list)',
                                       setup="from __main__ import brute_force, random_list", number=1)
    temp.append(bruteforce_runtime)
    print temp
    """
    """
    temp = []
    for y in [2, 3, 4, 5, 6, 7]:
        allowed_scores = [i for i in range(-10, 11)]
        random_list = [random.choice(allowed_scores) for x in range(10**y)]
        divide_runtime = timeit.timeit('divide_conquer(random_list)',
                                           setup="from __main__ import divide_conquer, random_list", number=1)
        temp.append(divide_runtime)
    print temp
    """
    """
    temp = []
    for y in [2, 3, 4, 5, 6, 7, 8]:
        allowed_scores = [i for i in range(-10, 11)]
        random_list = [random.choice(allowed_scores) for x in range(10**y)]
        linear_runtime = timeit.timeit('linear(random_list)',
                                       setup="from __main__ import linear, random_list", number=1)
        temp.append(linear_runtime)
    print temp
    """


