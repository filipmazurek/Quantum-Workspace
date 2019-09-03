import numpy as np
from grove.pyqaoa.maxcut_qaoa import maxcut_qaoa
from pyquil import get_qc
from pyquil.api import WavefunctionSimulator
from hidden_prints import hidden_prints
from rigetti_result_analysis import convert_result, plot_state_histogram

graph = [(0, 1), (0, 2), (0, 3)]

qvm = get_qc('4q-qvm')
with hidden_prints():
    maxcut_solver = maxcut_qaoa(graph=graph)
    betas, gammas = maxcut_solver.get_angles()

angles = np.hstack((betas, gammas))
param_prog = maxcut_solver.get_parameterized_program()
prog = param_prog(angles)  # gives the sequence of qubits and gates

measurements = qvm.run_and_measure(prog, trials=1000)  # simulate the program runs

counts = convert_result(measurements)
# plot_state_histogram(counts)

wavefunction = WavefunctionSimulator().wavefunction(prog)  # use wavefunction class to get theoretical measures
prob_dict = wavefunction.get_outcome_probs()  # extracts the probabilities of outcomes as a dict
print(prob_dict)
prob_dict.keys()  # these store the bitstring outcomes
