from numpy import kron as k
from numpy import subtract as sub
from numpy import add
from numpy import sqrt
from numpy import array
import numpy as np

pauli_Z = array([[1, 0], [0, -1]])
pauli_X = ([[0, 1], [1, 0]])
Id = ([[1, 0], [0, 1]])

# utility for 4 qubits
z0 = k(k(k(pauli_Z, Id), Id), Id)
z1 = k(k(k(Id, pauli_Z), Id), Id)
z2 = k(k(k(Id, Id), pauli_Z), Id)
z3 = k(k(k(Id, Id), Id), pauli_Z)
I = k(k(k(Id, Id), Id), Id)


def create_penalty_for_range_02():
    weight = -100

    z_term = z0 * -100
    all_ones_term = sub(I * -50, z0 * -50)

    z_term = z_term.dot(z2)
    all_ones_term = all_ones_term.dot(sub((I * .5), (z2 * .5)))

    cost_op = sub(sub(I * -100, z_term), all_ones_term)

    # print(cost_op)

    return cost_op


def create_penalty_for_range_13():
    weight = -100

    z_term = z1 * -100
    all_ones_term = sub(I * -50, z1 * -50)

    z_term = z_term.dot(z3)
    all_ones_term = all_ones_term.dot(sub((I * .5), (z3 * .5)))

    cost_op = sub(sub(I * -100, z_term), all_ones_term)

    # print(cost_op)

    return cost_op


def create_penalty_for_range_01():
    weight = -100

    z_term = z0 * -100
    all_ones_term = sub(I * -50, z0 * -50)

    z_term = z_term.dot(z1)
    all_ones_term = all_ones_term.dot(sub((I * .5), (z1 * .5)))

    cost_op = sub(sub(I * -100, z_term), all_ones_term)

    # print(cost_op)

    return cost_op


def create_penalty_for_range_23():
    weight = -100

    z_term = z2 * -100
    all_ones_term = sub(I * -50, z2 * -50)

    z_term = z_term.dot(z3)
    all_ones_term = all_ones_term.dot(sub((I * .5), (z3 * .5)))

    cost_op = sub(sub(I * -100, z_term), all_ones_term)

    # print(cost_op)

    return cost_op


def create_penalty_operators_for_qubit_range(self, range_of_qubits):
    from pyquil.paulis import PauliTerm, PauliSum
    cost_operators = []
    weight = -100 * np.max(self.distance_matrix)
    for i in range_of_qubits:
        if i == range_of_qubits[0]:
            z_term = PauliTerm("Z", i, weight)
            all_ones_term = PauliTerm("I", 0, 0.5 * weight) - PauliTerm("Z", i, 0.5 * weight)
        else:
            z_term = z_term * PauliTerm("Z", i)
            all_ones_term = all_ones_term * (PauliTerm("I", 0, 0.5) - PauliTerm("Z", i, 0.5))

    z_term = PauliSum([z_term])
    cost_operators.append(PauliTerm("I", 0, weight) - z_term - all_ones_term)

    return cost_operators


same_city_penalty = add(add(add(create_penalty_for_range_02(),
                                create_penalty_for_range_13()),
                            create_penalty_for_range_01()),
                        create_penalty_for_range_23())
# print(same_city_penalty)


def CRX(theta):
    CRX = array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, np.cos(theta / 2), -1j * np.sin(theta / 2)],
        [0, 0, -1j * np.sin(theta / 2), np.cos(theta / 2)]
        ])
    return CRX


Z = [[1, 0], [0, -1]]
X = [[0, 1], [1, 0]]
Id = [[1, 0], [0, 1]]
H = 1/sqrt(2) * array([[1, 1], [1, -1]])
S = array([[1, 0], [0, 1j]])

CNOT = array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]])

CNOT_rev = array([
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]])

qubits = array([[1], [2], [3], [4]])

# print(CNOT_rev @ qubits)

