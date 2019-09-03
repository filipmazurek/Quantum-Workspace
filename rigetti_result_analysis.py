"""
Filip Mazurek - 9/1/2019

Utility to make analyzing results from pyquil easier
"""

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt


def convert_result(measurements):
    """
    :param measurements: results from run_and_measure() using pyquil
    :return: Counter specifying frequency of results per trial
    """
    num_qubits = len(measurements)
    num_trials = len(measurements[0])

    results_per_trial = [[-1 for y in range(num_qubits)] for x in range(num_trials)]

    for i in range(len(measurements)):  # number of trials
        for j in range(len(measurements[0])):  # number of qubits
            results_per_trial[j][i] = measurements[i][j]

    # a hack so that we can use the Counter. Counter will take in tuples, but not lists
    tupled_result = [tuple(result) for result in results_per_trial]

    return Counter(tupled_result).most_common()


def plot_state_histogram(states_with_probs):
    states = np.array(states_with_probs)[:,0]
    probs = np.array(states_with_probs)[:,1].astype(float)
    n = len(states_with_probs)
    plt.barh(range(n), probs, tick_label=states)
    plt.show()