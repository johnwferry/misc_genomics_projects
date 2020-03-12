from compsci260lib import *
from UltrametricAdditive import is_ultrametric, is_additive, is_almost_equal


# Global variables available to all functions
# Get the aligned DNA sequences:
fasta_dict = get_fasta_dict('CoV.aligned.fasta')
seq_names = {1: 'AAK83356.1', 2: 'AAL57308.1', 3: 'AAL80031',
             4: 'AAP41037.1', 5: 'AAR01015.1', 6: 'AAV49723',
             7: 'ABD75513.1', 8: 'AFS88936.1', 9: 'ATW75480.1',
             10: 'AUM60024.1', 11: "CAA01637.1", 12: 'CAA71056.1',
             13: 'YP_209233.1'}


def build_tree():
    # Compute the pairwise distances:
    dist = {}
    for i in range(1, len(list(seq_names.keys()))):
        for j in range(i+1, len(list(seq_names.keys())) + 1):
            # Write the code to complete the compute_dist function
            dist[str(i) + "," + str(j)] = compute_dist(i, j)

    # uncomment the following line to print out the distance dictionary
    #print_dist(dist)

    # Check if the distance is ultrametric.
    if is_ultrametric(dist, threshold=0.0275):
        print("\nThe distance is ultrametric.")
    else:
        print("\nThe distance is not ultrametric.")

    # Check if the distance is additive.
    if is_additive(dist, threshold=0.0275):
        print("\nThe distance is additive.\n")
    else:
        print("\nThe distance is not additive.\n")

    # ------------- Implementation of the Neighbor Joining algorithm ----------
    # Keeping track of variable names:

    # dist              - dictionary containing the computed pair-wise
    #                     distances between nodes
    # node_list         - array to contain the list of current nodes, while
    #                     iterating to build the tree
    # avg_dist          - dictionary to contain the averaged distances for all
    #                     current nodes
    # D                 - dictionary to contain the adjusted distances between
    #                     nodes
    # newick            - dictionary to maintain the Newick representation of
    #                     each leaf node or subtree

    # ------------- Initialization: -------------------------------------------

    node_list = list(seq_names.keys())  # L = the list of current nodes
    newick = {}
    for i in range(1, len(node_list)+1):
        newick[i] = "" + str(i)

    avg_dist = {}  # averaged distances (the r terms) for all current nodes

    for i in range(1, len(node_list) + 1):
        avg_dist[i] = 0
        for j in range(1, len(node_list) + 1):
            if i != j:
                # if i < j, then  avg_dist[i] = dist["i,j"]
                # else avg_dist[i] = dist["j,i"]
                avg_dist[i] += dist["%d,%d" %
                                    (i, j)] if i < j else dist["%d,%d" %
                                                               (j, i)]
        avg_dist[i] = avg_dist[i] / (len(node_list) - 2)

    max_node = len(node_list)  # the maximum key used for a node

    # -------------- Iteration to build the tree --------------------

    # As long as the list contains more than 2 nodes, iterate
    while len(node_list) > 2:

        # ---------- Begin your code -------------

        # Compute the 'adjusted' distances between nodes using the original
        # distances (from dist)
        # and averaged distances (from avg_dist)

        # Let D be the dict to contain 'adjusted' distances
        D = {}

        # Loop through each pair of nodes as entered in the dist dict
        # Use the entries from the avg_dist and calculate entries for the D
        # dict


        # Pick the pair i,j in node_list for which adjusted distance D_ij is
        # minimal.
        # You may find the function two_min_in_dict helpful.
        (i, j) = None, None  # Replace with your pair

        # Define a new node k and set dist[m,k] for all nodes m in node_list
        # as (dist[i,m] + dist[j,m] - dist[i,j]) / 2

        # max_node had been earlier set to the largest key used for a node
        k = max_node + 1
        max_node += 1
        m = 0
        for ind in range(len(node_list)):
            m = node_list[ind]

        # ---------- End your code -------------

        # Add the new node k to the Newick format representation of the tree
        # with edges of lengths
        #dik = (dist[i,j] + avg_dist[i] - avg_dist[j])/2
        # djk = dist[i,j]-d[i,k],
        # joining k to i and j

        d_ik = (dist["%d,%d" % (i, j)] + avg_dist[i] - avg_dist[j]) / 2
        d_jk = dist["%d,%d" % (i, j)] - d_ik
        newick[k] = "(" + newick[i] + ":" + "%.4f" % (d_ik) + \
            "," + newick[j] + ":" + "%.4f" % (d_jk) + ")"

        # Remove i and j from node_list and add k
        temp = []
        for idx in range(len(node_list)):
            if node_list[idx] != i and node_list[idx] != j:
                temp.append(node_list[idx])
        temp.append(k)

        node_list = list(temp)

        # Update the r terms
        if len(node_list) > 2:
            avg_dist[k] = 0
            for ind in range(len(node_list) - 1):
                m = node_list[ind]
                avg_dist[m] = avg_dist[m] * (len(node_list) - 1)
                avg_dist[m] -= dist["%d,%d" %
                                    (m, i)] \
                    if m < i else dist["%d,%d" % (i, m)]
                avg_dist[m] -= dist["%d,%d" %
                                    (m, j)] \
                    if m < j else dist["%d,%d" % (j, m)]
                avg_dist[m] += dist["%d,%d" % (m, k)]

                avg_dist[m] /= (len(node_list) - 2)
                avg_dist[k] += dist["%d,%d" % (m, k)]

            avg_dist[k] = avg_dist[k] / (len(node_list) - 2)

        # Remove any elements from the dict that contain nodes i or j
        delete_from_dict(dist, i, j)
        delete_from_dict(D, i, j)
        delete_from_dict(newick, i, j)
        delete_from_dict(avg_dist, i, j)

    print("The Newick representation for the phylogenetic tree is:\n(" +
          newick[node_list[0]] + ":" + "%.4f" % (list(dist.values())[0]) +
          "," + newick[node_list[1]] + ":0);\n")


def compute_dist(i, j):
    """Returns the distance between two sequences. The distance is computed
    as the ratio between the number of mismatches and the total number
    of matches or mismatches."""

    seq1 = fasta_dict[seq_names[i]]
    seq2 = fasta_dict[seq_names[j]]
    mismatch = 0
    matchORmismatch = 0

    for i in range(0, len(seq1)):
        # ---------- Begin your code -------------

        # Fill in your code to count the number of matches and mismatches.
        # Ignore gaps in either sequence.

        if i >= len(seq2)  or seq1[i] != seq2[i]:
            mismatch +=1
        matchORmismatch +=1



        # ---------- End your code -------------
        pass

    return float(mismatch)/matchORmismatch


def delete_from_dict(dictionary, i, j):
    """Deletes the dict entries with keys that contain i or j.  Doesn't return
    anything; just changes the dictionary."""
    for k in list(dictionary.keys()):
        ks = [int(_) for _ in str(k).split(',')]
        if i in ks or j in ks:
            del dictionary[k]


def min_in_dict(dictionary):
    """Returns the key associated with the minimum value in the dict."""
    import operator
    return min(iter(dictionary.items()), key=operator.itemgetter(1))[0]


def two_min_in_dict(dictionary):
    import operator
    sorted_dict = sorted(iter(dictionary.items()), key=operator.itemgetter(1))
    element = sorted_dict[0][0]  # get the first element of the tuple
    (i, j) = element.split(",")
    return (int(i), int(j))


def print_dist(dist):
    '''print the distance store in the given dict structure'''
    idx = [int(i.split(',')[0]) for i in list(dist.keys())]
    idx.extend([int(i.split(',')[1]) for i in list(dist.keys())])
    max_idx = max(idx)

    print(' \t', end=' ')
    for col in range(2, max_idx + 1):
        print('%d\t' % col, end=' ')

    for row in range(1, max_idx):
        print('%d\t' % row, end=' ')
        for col in range(2, row+1):
            print(' \t', end=' ')
        for col in range(row+1, max_idx + 1):
            print('%.4f\t' % dist[str(row)+','+str(col)], end=' ')


if __name__ == '__main__':
    build_tree()
