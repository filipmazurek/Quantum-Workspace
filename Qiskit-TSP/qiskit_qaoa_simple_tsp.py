"""
Filip Mazurek 9/10/2019
A simple version of the travelling salesman. This is only to find any Hamiltonian path

Made with help from the qiskit community tutorial on qaoa: qml_mooc
https://github.com/Qiskit/qiskit-community-tutorials/blob/41dfb73df77595c5d5078164ddecfd915e24897b/awards/teach_me_quantum_2018/qml_mooc/07_Variational%20Circuits.ipynb
"""
from qiskit.quantum_info import Pauli
from qiskit.aqua import Operator
import numpy as np
from functools import partial, reduce
from qiskit.aqua.components.initial_states import Custom
from qiskit import QuantumRegister, Aer, BasicAer, execute
from scipy.optimize import minimize
from qiskit.aqua.translators.ising import max_cut, tsp

np.set_printoptions(precision=3, suppress=True)


def pauli_x(qubit, coeff, n_q):
    eye = np.eye(n_q)
    return Operator([[coeff, Pauli(np.zeros(n_q), eye[qubit])]])


def pauli_z(qubit, coeff, n_q):
    eye = np.eye(n_q)
    return Operator([[coeff, Pauli(eye[qubit], np.zeros(n_q))]])


def product_pauli_z(q1, q2, coeff, n_q):
    eye = np.eye(n_q)
    return Operator([[coeff, Pauli(eye[q1], np.zeros(n_q)) * Pauli(eye[q2], np.zeros(n_q))]])


def evolve(hamiltonian, angle, quantum_registers):
    return hamiltonian.evolve(None, angle, 'circuit', 1,
                              quantum_registers=quantum_registers,
                              expansion_mode='suzuki',
                              expansion_order=3)


def create_circuit(qr, gamma, beta, p, m_H, c_H, init_circ):
    circuit_evolv = reduce(lambda x, y: x + y, [evolve(m_H, beta[i], qr) + evolve(c_H, gamma[i], qr)
                                                for i in range(p)])
    circuit = init_circ + circuit_evolv
    return circuit


def neg_evaluate_circuit(gamma_beta, qr, p, m_H, c_H, init_circ):
    n = len(gamma_beta)//2
    circuit = create_circuit(qr, gamma_beta[:n], gamma_beta[n:], p, m_H=m_H,  c_H=c_H, init_circ=init_circ)
    return -1 * np.real(c_H.eval("matrix", circuit, Aer.get_backend('statevector_simulator'))[0])


def create_penalty_operators_for_qubit_range(range_of_qubits, dist_mat):
    penalty_weight = 100 * np.max(dist_mat)



    pass

# def create_penalty_operators_for_qubit_range(self, range_of_qubits):
#     cost_operators = []
#     weight = -100 * np.max(self.distance_matrix)
#     for i in range_of_qubits:
#         if i == range_of_qubits[0]:
#             z_term = PauliTerm("Z", i, weight)
#             all_ones_term = PauliTerm("I", 0, 0.5 * weight) - PauliTerm("Z", i, 0.5 * weight)
#         else:
#             z_term = z_term * PauliTerm("Z", i)
#             all_ones_term = all_ones_term * (PauliTerm("I", 0, 0.5) - PauliTerm("Z", i, 0.5))
#
#     z_term = PauliSum([z_term])
#     cost_operators.append(PauliTerm("I", 0, weight) - z_term - all_ones_term)
#
#     return cost_operators

def main():
    # graph of city coordinates
    cities = np.array([[0, 0], [0, 1]])
    num_cities = len(cities)
    num_qubits = num_cities ** 2

    # algorithm properties
    p = 2  # number of time steps
    beta = np.random.uniform(0, np.pi * 2, p)
    gamma = np.random.uniform(0, np.pi * 2, p)

    # create matrix of distances between cities
    distance_mat = tsp.calc_distance(cities).w  # note that this method does integer distances

    # create mixing Hamiltonian. A city may or may not be visited in a timestep
    mixing_hamiltonian = reduce(lambda x, y: x + y,
                                [pauli_x(i, 1, num_qubits) for i in range(num_qubits)])

    for t in range(num_cities):
        range_of_qubits = list(range(t * num_cities, (t + 1) * num_cities))
        print(range_of_qubits)

        for i in range_of_qubits:
            print(i)
    pass



main()
