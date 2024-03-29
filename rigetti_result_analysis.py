"""
Filip Mazurek - 9/1/2019

Utility to make analyzing results from pyquil easier
"""

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt


def convert_result(measurements):
    # TODO: please note the endian-ness and how that affects everything else. Right now least significant bit is on left
    # TODO: rewrite using numpy.vstack instead of this mess
    """
    :param measurements: results from run_and_measure() using pyquil
    :return: Counter object. May need to use most_common function.
    """
    num_qubits = len(measurements)
    num_trials = len(measurements[0])

    results_per_trial = [[-1 for y in range(num_qubits)] for x in range(num_trials)]

    for i in range(len(measurements)):  # number of trials
        for j in range(len(measurements[0])):  # number of qubits
            results_per_trial[j][i] = measurements[i][j]

    # a hack so that we can use the Counter. Counter will take in tuples, but not lists
    tupled_result = [tuple(result) for result in results_per_trial]

    return Counter(tupled_result)


def plot_state_histogram(states_with_probs):
    states = np.array(states_with_probs)[:,0]
    probs = np.array(states_with_probs)[:,1].astype(float)
    n = len(states_with_probs)
    plt.barh(range(n), probs, tick_label=states)
    plt.show()


def error_binary_state_to_points_order(binary_state):
    """
    A modification on MS's original function. This will sort out the erroneous results as such, and will
    only keep the results which make sense
    Transforms the the order of points from the binary representation: [1,0,0,0,1,0,0,0,1],
    to the standard one: [0, 1, 2]

    Transforms [1,1,0,0] to erroneous

    NOTE: This method assumes that the binary state is a matrix row-by-row.
    :param binary_state:
    :return: standard lists
    """
    points_order = []
    number_of_points = int(np.sqrt(len(binary_state)))
    column_points = []
    error_rep = [-1]
    for p in range(number_of_points):
        row_done = False
        for j in range(number_of_points):
            if binary_state[number_of_points * p + j] == 1:
                if row_done:  # there is already a 1 in this row
                    return error_rep
                elif p in column_points:  # there is already a 1 in this column
                    return error_rep
                else:
                    points_order.append(j)
                    row_done = True
                    column_points.append(p)

    if len(points_order) != number_of_points:  # there were not enough ones
        return error_rep
    else:
        return points_order


def tsp_convert_raw_to_order(sampling_results):  # TODO: check for usage and delete
    """
    :param raw_sampling: the result of the quantum computer running the tsp algorithm
    :return: show which sensible results are left. Discard nonsensical answers (two cities at the same time, etc.)
    """
    all_solutions = sampling_results.keys()
    naive_distribution = {}
    for sol in all_solutions:
        points_order_solution = error_binary_state_to_points_order(sol)
        if tuple(points_order_solution) in naive_distribution.keys():  # Can this ever be true?
            naive_distribution[tuple(points_order_solution)] += sampling_results[sol]
        else:
            naive_distribution[tuple(points_order_solution)] = sampling_results[sol]

    pass
