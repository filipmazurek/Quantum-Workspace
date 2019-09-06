"""
Tutorial on variational circuits. This completely covers qaoa. Just replace the Hamiltonians for ones that correspond
to the problem being solved. May also need to make more Pauli operators to be useful.
https://github.com/Qiskit/qiskit-community-tutorials/blob/41dfb73df77595c5d5078164ddecfd915e24897b/awards/teach_me_quantum_2018/qml_mooc/07_Variational%20Circuits.ipynb
"""
import itertools
import numpy as np
from functools import partial, reduce
from qiskit import Aer, BasicAer, QuantumRegister, execute
from qiskit.quantum_info import Pauli
from qiskit.aqua import Operator
from qiskit.aqua.components.initial_states import Custom
from scipy.optimize import minimize
np.set_printoptions(precision=3, suppress=True)

N_QUBITS = 2


def pauli_x(qubit, coeff):
    eye = np.eye(N_QUBITS)
    return Operator([[coeff, Pauli(np.zeros(N_QUBITS), eye[qubit])]])


def pauli_z(qubit, coeff):
    eye = np.eye(N_QUBITS)
    return Operator([[coeff, Pauli(eye[qubit], np.zeros(N_QUBITS))]])


def product_pauli_z(q1, q2, coeff):
    eye = np.eye(N_QUBITS)
    return Operator([[coeff, Pauli(eye[q1], np.zeros(N_QUBITS)) * Pauli(eye[q2], np.zeros(N_QUBITS))]])


def evolve(hamiltonian, angle, quantum_registers):
    return hamiltonian.evolve(None, angle, 'circuit', 1,
                              quantum_registers=quantum_registers,
                              expansion_mode='suzuki',
                              expansion_order=3)


def create_circuit(qr, gamma, beta, p):
    circuit_evolv = reduce(lambda x,y: x+y, [evolve(Hm, beta[i], qr) + evolve(Hc, gamma[i], qr)
                                             for i in range(p)])
    circuit = circuit_init + circuit_evolv
    return circuit


def evaluate_circuit(gamma_beta, qr, p):
    n = len(gamma_beta)//2
    circuit = create_circuit(qr, gamma_beta[:n], gamma_beta[n:], p)
    return np.real(Hc.eval("matrix", circuit, Aer.get_backend('statevector_simulator'))[0])


# Create the driving (mixing) Hamiltonian
Hm = reduce(lambda x, y: x+y,
            [pauli_x(i, 1) for i in range(N_QUBITS)])
Hm.to_matrix()

# Create the cost Hamiltonian
J = np.array([[0, 1], [0, 0]])
Hc = reduce(lambda x, y: x + y,
            [product_pauli_z(i, j, -J[i, j])
             for i, j in itertools.product(range(N_QUBITS), repeat=2)])
Hc.to_matrix()

n_iter = 10  # number of iterations of the optimization procedure
p = 2  # number of operators
beta = np.random.uniform(0, np.pi*2, p)  # set up random angles to start
gamma = np.random.uniform(0, np.pi*2, p)

# Initializing the initial state to an equal superposition of all states. "1" is given as the possibility for all states
init_state_vect = [1 for i in range(2**N_QUBITS)]
init_state = Custom(N_QUBITS, state_vector=init_state_vect)

# Prepare the actual initial state with a quantum register
qr = QuantumRegister(N_QUBITS)
circuit_init = init_state.construct_circuit('circuit', qr)

# Create as a partial so that when minimizing we don't change the values of the qr or of p.
evaluate = partial(evaluate_circuit, qr=qr, p=p)

# Minimize evaluate_circuit over gamma and beta only.
# We are minimizing so that we can have the minimum cost Hamiltonian
# This is the "variational" part of the algorithm. We are running the quantum circuit multiple times to find the minimum
result = minimize(evaluate, np.concatenate([gamma, beta]), method='L-BFGS-B')
print(result)

# The above found the minimum beta and gamma. So we run the quantum circuit with these angles

circuit = create_circuit(qr, result['x'][:p], result['x'][p:], p)

backend = BasicAer.get_backend('statevector_simulator')
job = execute(circuit, backend)
state = np.asarray(job.result().get_statevector(circuit))
print(np.absolute(state))
print(np.angle(state))

Z0 = pauli_z(0, 1)
Z1 = pauli_z(1, 1)

print(Z0.eval("matrix", circuit, Aer.get_backend('statevector_simulator'))[0])
print(Z1.eval("matrix", circuit, Aer.get_backend('statevector_simulator'))[0])
