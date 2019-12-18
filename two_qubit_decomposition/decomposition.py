import numpy as np
from math import isclose
from cirq.linalg.decompositions import kron_factor_4x4_to_2x2s, so4_to_magic_su2s
from cirq.linalg.predicates import is_special_orthogonal


def factor_kronecker(mat: np.array) -> [np.array]:
    pass


pauli_Z = np.array([[1, 0], [0, -1]])
pauli_X = np.array([[0, 1], [1, 0]])
I2 = np.eye(2)
I4 = np.eye(4)
H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
S = np.array([[1, 0], [0, 1j]])
CNOT1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
CNOT2 = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
magic_basis = CNOT2 @ np.kron(I2, H) @ np.kron(I2, S) @ np.kron(S, I2)

# U that are real and det = 1
# U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ CNOT1 @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2) @ CNOT2
# U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2)
# U = np.kron(pauli_X, I2)

# U that are real and det = -1
U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ CNOT1 @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2)


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
            print("Nonreal matrix")

# print("U\n")
# print(U)

print(is_special_orthogonal(U))

print("Determinant = ", np.around(np.linalg.det(U), 2))

mUm = magic_basis @ U @ magic_basis.conjugate().transpose() @ np.kron(I2, pauli_Z) @ CNOT1 @ CNOT2 @ CNOT1

tuple = kron_factor_4x4_to_2x2s(mUm)

# print(tuple[1])
# print(tuple[2])

print()
print(U)
# print()
# print(np.around(magic_basis.conjugate().transpose() @ U @ magic_basis, 2))
print()
print(np.around(magic_basis.conjugate().transpose() @ np.kron(tuple[1], tuple[2]) * tuple[0]  @ CNOT1 @ CNOT2 @ CNOT1 @ np.kron(I2, pauli_Z) @ magic_basis, 2))


# tuple = so4_to_magic_su2s(U)
# Mag = np.array([[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]])
# print(Mag.conjugate().transpose() @ np.kron(tuple[0], tuple[1]) @ Mag)

print(np.kron(pauli_Z, pauli_Z))
print(np.kron(pauli_X, pauli_X))
