"""
Filip Mazurek 9/10/2019
A simple version of the travelling salesman. This is only to find any Hamiltonian path

Made with help from the qiskit community tutorial on qaoa: qml_mooc
https://github.com/Qiskit/qiskit-community-tutorials/blob/
41dfb73df77595c5d5078164ddecfd915e24897b/awards/teach_me_quantum_2018/qml_mooc/07_Variational%20Circuits.ipynb
"""
from qiskit.quantum_info import Pauli
from qiskit.aqua import Operator
import numpy as np
from functools import partial, reduce
from qiskit.aqua.components.initial_states import Custom
from qiskit import QuantumRegister, Aer, BasicAer, execute
from scipy.optimize import minimize
from qiskit.aqua.translators.ising import max_cut, tsp
from results_visualization import list_to_easier_vis
np.set_printoptions(precision=3, suppress=True)


def pauli_i(coeff, n_q):
    id_pauli = Pauli(np.zeros(n_q), np.zeros(n_q))
    return Operator([[coeff, id_pauli]])


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
    return np.real(c_H.eval("matrix", circuit, Aer.get_backend('statevector_simulator'))[0])


def create_weights_cost_operators(num_cities, num_qubits, dist_mat):
    cost_operator = None

    for i in range(num_cities):
        for j in range(i, num_cities):
            for t in range(num_cities - 1):
                weight = dist_mat[i][j] / 2
                if dist_mat[i][j] != 0:
                    qubit_1 = t * num_cities + i
                    qubit_2 = (t + 1) * num_cities + j
                    if cost_operator is None:
                        cost_operator = pauli_i(weight, num_qubits) - \
                                        product_pauli_z(qubit_1, qubit_2, weight, num_qubits)
                    else:
                        cost_operator += pauli_i(weight, num_qubits) - \
                                         product_pauli_z(qubit_1, qubit_2, weight, num_qubits)
    return cost_operator


def create_penalty_operators_for_bilocation(num_cities, distance_mat, num_qubits):
    # TODO: big problems here. It likes position 1010 WAAY too much (= 0.88) (in two city case)
    penalty_operators = None
    for t in range(num_cities):  # adding penalty for being in multiple cities at the same time point
        range_of_qubits = list(range(t * num_cities, (t + 1) * num_cities))
        print(range_of_qubits)
        if penalty_operators is None:
            penalty_operators = create_penalty_operators_for_qubit_range(range_of_qubits, distance_mat, num_qubits)
        else:
            penalty_operators += create_penalty_operators_for_qubit_range(range_of_qubits, distance_mat, num_qubits)

    return penalty_operators


def create_penalty_operators_for_repetition(num_cities, distance_mat, num_qubits):
    # TODO: big problems here. It likes position 1100 WAAY too much (= 0.88) (in two city case)
    penalty_operators = None
    for i in range(num_cities):  # add penalty for visiting the same city more than once
        range_of_qubits = list(range(i, num_cities ** 2, num_cities))
        print(range_of_qubits)
        if penalty_operators is None:
            penalty_operators = create_penalty_operators_for_qubit_range(range_of_qubits, distance_mat, num_qubits)
        else:
            penalty_operators += create_penalty_operators_for_qubit_range(range_of_qubits, distance_mat, num_qubits)
    return penalty_operators


def create_penalty_operators_for_qubit_range(range_of_qubits, dist_mat, n_q):
    penalty_weight = 100 * np.max(dist_mat)
    cost_operators = None
    for i in range_of_qubits:
        if i == range_of_qubits[0]:
            z_term = pauli_z(qubit=i, coeff=penalty_weight, n_q=n_q)
            all_ones_term = pauli_i(coeff=.5 * penalty_weight, n_q=n_q) - pauli_z(qubit=i, coeff=0.5 * penalty_weight, n_q=n_q)
        else:
            z_term = z_term * pauli_z(qubit=i, coeff=1, n_q=n_q)
            all_ones_term = all_ones_term * (pauli_i(coeff=.5, n_q=n_q) - pauli_z(qubit=i, coeff=0.5, n_q=n_q))

        if cost_operators is None:
            cost_operators = pauli_i(penalty_weight, n_q) - z_term - all_ones_term
        else:
            cost_operators += pauli_i(penalty_weight, n_q) - z_term - all_ones_term

    return cost_operators


def main(run_mode):
    # graph of city coordinates
    cities = np.array([[0, 0], [0, 1]])  # coordinates of the cities
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

    # penalty_operators = create_weights_cost_operators(num_cities=num_cities, num_qubits=num_qubits,
    #                                                   dist_mat=distance_mat)
    # penalty_operators += create_penalty_operators_for_bilocation(num_qubits=num_qubits, num_cities=num_cities,
    #                                                              distance_mat=distance_mat)
    penalty_operators = create_penalty_operators_for_repetition(num_qubits=num_qubits, num_cities=num_cities,
                                                                 distance_mat=distance_mat)

    print(penalty_operators)
    cost_hamiltonian = penalty_operators

    # circuit initial state vector. All states in equal superposition
    init_state_vect = [1 for i in range(2 ** num_qubits)]
    init_state = Custom(num_qubits, state_vector=init_state_vect)

    # initialize quantum circuit
    qr = QuantumRegister(num_qubits)
    init_circ = init_state.construct_circuit('circuit', qr)

    # find optimal beta and gamma
    evaluate = partial(neg_evaluate_circuit, qr=qr, p=p, m_H=mixing_hamiltonian, c_H=cost_hamiltonian,
                       init_circ=init_circ)
    print("Looking for optimal beta and gamma")
    # TODO: maybe we should use a different or faster method of finding the min? Super long even with two cities
    result = minimize(evaluate, np.concatenate([gamma, beta]), method='L-BFGS-B')
    # result = minimize(evaluate, np.concatenate([gamma, beta]))

    print(result)

    # now use the result of the gathered angles to find the answer
    circuit = create_circuit(qr, result['x'][:p], result['x'][p:], p, m_H=mixing_hamiltonian, c_H=cost_hamiltonian,
                             init_circ=init_circ)

    if run_mode == "IBM quantum":
        import secrets
        from qiskit import IBMQ
        from qiskit.providers.ibmq import least_busy

        provider = IBMQ.enable_account(secrets.IBM_TOKEN)
        large_enough_devices = provider.backends(filters=lambda x: x.configuration().n_qubits > 4 and
                                                                   not x.configuration().simulator)
        backend = least_busy(large_enough_devices)
        print("This will be running on the IBM device " + backend.name())

    else:
        print("Preparing to run on local simulator")
        backend = BasicAer.get_backend('statevector_simulator')

    job = execute(circuit, backend)
    state = np.asarray(job.result().get_statevector(circuit))
    print(list_to_easier_vis(np.absolute(state)))


main(run_mode="sim")
