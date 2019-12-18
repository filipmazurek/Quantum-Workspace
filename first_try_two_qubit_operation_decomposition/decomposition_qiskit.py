from qiskit import QuantumCircuit, Aer, execute
import numpy as np
from math import isclose
from scipy.optimize import minimize
import cirq.optimizers


pauli_Z = np.array([[1, 0], [0, -1]])
pauli_X = np.array([[0, 1], [1, 0]])
I2 = np.eye(2)
I4 = np.eye(4)
H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
S = np.array([[1, 0], [0, 1j]])
CNOT1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
CNOT2 = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
magic_gate = CNOT2 @ np.kron(I2, H) @ np.kron(I2, S) @ np.kron(S, I2)  # magic gate verified as correct


def round_matrix(mat):
    rounded_result = np.zeros((4, 4))

    for i in range(4):
        for j in range(4):
            rounded_result[i][j] = round(mat[i][j].real, 1) + round(mat[i][j].imag, 2) * 1j
    return rounded_result


def SO4_circuit(params: list):
    # Construct an empty quantum circuit
    circ = QuantumCircuit(2)

    circ.rz(np.pi/2, 0)
    circ.rz(np.pi/2, 1)
    circ.ry(np.pi/2, 1)

    # Either the above section or the below. Only differ by a cancelled phase
    # circ.s(0)
    # circ.s(1)
    # circ.h(1)

    circ.cx(1, 0)

    circ.rz(params[0], 0)
    circ.ry(params[1], 0)
    circ.rz(params[2], 0)

    circ.rz(params[3], 1)
    circ.ry(params[4], 1)
    circ.rz(params[5], 1)

    circ.cx(1, 0)

    circ.ry(-np.pi/2, 1)
    circ.rz(-np.pi/2, 0)
    circ.rz(-np.pi/2, 1)

    # Either the above section or the below. Only differ by a cancelled phase
    # circ.h(1)
    # circ.sdg(1)
    # circ.sdg(0)

    # print(circ)
    # Select the UnitarySimulator from the Aer provider
    simulator = Aer.get_backend('unitary_simulator')

    # Execute and get counts
    result = execute(circ, simulator).result()
    unitary = result.get_unitary(circ)
    # print("Circuit unitary:\n", unitary)
    return unitary


def cost_function_SO4(params: list):
    cost = 0
    SO4 = SO4_circuit(params)

    identity_goal = SO4 @ np.linalg.inv(U)
    for i in range(4):
        for j in range(4):
            cost += abs(identity_goal[i][j] - I4[i][j])

    return cost


def SU4_circuit(params: list):
    circ = QuantumCircuit(2)

    circ.rz(params[0], 0)
    circ.ry(params[1], 0)
    circ.rz(params[2], 0)

    circ.rz(params[3], 1)
    circ.ry(params[4], 1)
    circ.rz(params[5], 1)

    circ.cx(1, 0)

    circ.rz(params[6], 0)
    circ.ry(params[7], 1)

    circ.cx(0, 1)

    circ.ry(params[8], 1)

    circ.cx(1, 0)

    circ.rz(params[9], 0)
    circ.ry(params[10], 0)
    circ.rz(params[11], 0)

    circ.rz(params[12], 1)
    circ.ry(params[13], 1)
    circ.rz(params[14], 1)

    simulator = Aer.get_backend('unitary_simulator')

    # Execute and get counts
    result = execute(circ, simulator).result()
    unitary = result.get_unitary(circ)
    # print("Circuit unitary:\n", unitary)
    return unitary

def cost_function_SU4(params: list):
    cost = 0
    SU4 = SU4_circuit(params)

    identity_goal = SU4 @ np.linalg.inv(U)
    for i in range(4):
        for j in range(4):
            cost += abs(identity_goal[i][j] - I4[i][j])

    return cost


U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ CNOT1 @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2) @ CNOT2
# U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2)
# U = CNOT1 @ CNOT2 @ CNOT1
# U = np.kron(pauli_X, I2)

print(U)
print(np.linalg.det(U))

# given an input U
# check that it is a valid 2 qubit operation
try:
    np.linalg.inv(U)
except np.linalg.LinAlgError:
    print("Not a valid unitary operation")

# check if it is all real. This lets us use 2 CNOTS
real = True
for row in U:
    for element in row:
        if np.imag(element) != 0:
            real = False
            print("Got a nonreal matrix")


optimized_result = None
circuit_result = None
# if all real, adjust method based on determinant (1 or -1)
if real:
    if isclose(np.linalg.det(U), 1):  # allows to use method with 2 CNOTS
        print()
        print("Real and det = 1")
        print()
        # parameters = np.array([np.pi - .5] * 6)
        parameters = np.array([np.pi - .5] * 15)
        # optimized_result = minimize(cost_function_SO4, parameters, method="L-BFGS-B")

        optimized_result = minimize(cost_function_SU4, parameters, method="L-BFGS-B")

        print("OPTIMIZED RESULT")
        print(optimized_result)
        print()

        # circuit_result = SO4_circuit(optimized_result['x'])
        circuit_result = SU4_circuit(optimized_result['x'])

print("CIRCUIT RESULT")
print(circuit_result)
print()

rounded_circuit_result = round_matrix(circuit_result)
print("ROUNDED CIRCUIT RESULT")
print(rounded_circuit_result)

# print(circuit_result)
print()
print("U")
print(U)

print()

product = circuit_result @ np.linalg.inv(U)
rounded_product = round_matrix(product)
print("ROUNDED PRODUCT")
print(rounded_product)


