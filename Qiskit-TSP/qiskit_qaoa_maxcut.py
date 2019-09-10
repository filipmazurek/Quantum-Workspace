"""
Filip Mazurek 9/6/2019
A custom qaoa solver in qiskit to solve the Max Cut problem

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


def main():
    # graph of edges to solve. If changing, update num_nodes as well
    graph = [(3, 0), (3, 1), (3, 2)]
    num_nodes = 4
    num_qubits = num_nodes

    # algorithm properties
    p = 2  # number of time steps
    beta = np.random.uniform(0, np.pi*2, p)
    gamma = np.random.uniform(0, np.pi*2, p)
    # n_iter = 10  # number of iterations of the optimization procedure

    # mixing hamiltonian. May pauli X each qubit (change their color for maxcut)
    mixing_hamiltonian = reduce(lambda x, y: x + y,
                                [pauli_x(i, 1, num_qubits) for i in range(num_qubits)])

    # identity operation. As a utility to get the correct cost hamiltonian
    id_pauli = Pauli(np.zeros(num_qubits), np.zeros(num_qubits))
    id_operation = Operator([[1, id_pauli]])

    # cost Hamiltonian: summation of pauli z products. To maximize the number of adjacent different nodes
    cost_hamiltonian = reduce(lambda x, y: x + y,
                              [id_operation - product_pauli_z(i, j, 1, n_q=num_qubits)
                               for i, j in graph])

    # circuit initial state vector. All states in equal superposition
    init_state_vect = [1 for i in range(2**num_qubits)]
    init_state = Custom(num_qubits, state_vector=init_state_vect)

    # initialize quantum circuit
    qr = QuantumRegister(num_qubits)
    init_circ = init_state.construct_circuit('circuit', qr)

    # find optimal beta and gamma
    evaluate = partial(neg_evaluate_circuit, qr=qr, p=p, m_H=mixing_hamiltonian, c_H=cost_hamiltonian, init_circ=init_circ)
    result = minimize(evaluate, np.concatenate([gamma, beta]), method='L-BFGS-B')
    print(result)

    # now use the result of the gathered angles to find the answer
    circuit = create_circuit(qr, result['x'][:p], result['x'][p:], p, m_H=mixing_hamiltonian, c_H=cost_hamiltonian, init_circ=init_circ)

    backend = BasicAer.get_backend('statevector_simulator')
    job = execute(circuit, backend)
    state = np.asarray(job.result().get_statevector(circuit))
    print(np.absolute(state))
    print(np.angle(state))


main()
