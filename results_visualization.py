from math import ceil


def list_to_easier_vis(prob_list):
    """
    Given a list of probabilities - where the position in the list corresponds to the qubit bitwise result.
    e.g. given [0.5, 0.5, .7, 0], return {00:.5, 01:.5, 10:.7, 11:0}
    :param list:
    :return:
    """
    prob_len = len(prob_list)
    prob_dict = {}
    for i in range(prob_len):
        prob_dict[binary(i, ceil(prob_len**.5))] = prob_list[i]

    return prob_dict


def binary(n, length):
    s = bin(n)
    # removing "0b" prefix
    s1 = s[2:]
    # add leading 0's
    while len(s1) < length:
        s1 = "0" + s1
    return s1
