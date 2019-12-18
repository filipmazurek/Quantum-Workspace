import numpy as np
from math import isclose
from scipy.optimize import minimize
import qiskit

# matrix = R_y(np.pi/2) @ np.linalg.inv(R_y(np.pi/2))
# for i in range(len(matrix)):
#     for j in range(len(matrix[i])):
#         matrix[i, j] = np.round(matrix[i][j], 14)

pauli_Z = np.array([[1, 0], [0, -1]])
pauli_X = np.array([[0, 1], [1, 0]])
I2 = np.eye(2)
I4 = np.eye(4)
H = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
S = np.array([[1, 0], [0, 1j]])
CNOT1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
CNOT2 = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])
magic_gate = CNOT2 @ np.kron(I2, H) @ np.kron(I2, S) @ np.kron(S, I2)  # magic gate verified as correct


def round_circuit(circ):
    rounded_circuit_result = np.zeros((4, 4))

    for i in range(4):
        for j in range(4):
            rounded_circuit_result[i][j] = round(circ[i][j].real, 1) + round(circ[i][j].imag, 2) * 1j
    return rounded_circuit_result


def R_y(theta):
    return np.array([[np.cos(theta/2), np.sin(theta/2)], [-np.sin(theta/2), np.cos(theta/2)]])


def R_z(alpha):
    return np.array([[np.exp(1j * alpha/2), 0], [0, np.exp(-1j * alpha/2)]])


S1 = R_z(np.pi/2)
S1_inv = np.linalg.inv(S1)
R1 = R_y(np.pi/2)
R1_inv = np.linalg.inv(R1)


def SO4_circuit(a_alpha, a_theta, a_beta, b_alpha, b_theta, b_beta):
    """
    Create the SO4 circuit to use to minimize
    :return:
    """
    # return np.kron(S1_inv, I2) @ np.kron(I2, S1_inv) @ np.kron(I2, R1_inv) @ CNOT2 \
    #        @ np.kron(I2, R_z(b_beta)) @ np.kron(I2, R_y(b_theta)) @ np.kron(I2, R_z(b_alpha)) \
    #        @ np.kron(R_z(a_beta), I2) @ np.kron(R_y(a_theta), I2) @ np.kron(R_z(a_alpha), I2) \
    #        @ CNOT2 @ np.kron(I2, R1) @ np.kron(I2, S1) @ np.kron(S1, I2)

    return np.linalg.inv(magic_gate) \
           @ np.kron(I2, R_z(b_beta)) @ np.kron(I2, R_y(b_theta)) @ np.kron(I2, R_z(b_alpha)) \
           @ np.kron(R_z(a_beta), I2) @ np.kron(R_y(a_theta), I2) @ np.kron(R_z(a_alpha), I2) \
           @ magic_gate


def cost_function_SO4(params: list):
    """
    Make a cost function to find all parameters
    :return:
    """
    cost = 0
    SO4 = SO4_circuit(params[0], params[1], params[2], params[3], params[4], params[5])

    for i in range(4):
        for j in range(4):
            cost += abs(SO4[i][j] - U[i][j])

    # identity_goal = SO4 @ np.linalg.inv(U)
    # for i in range(4):
    #     for j in range(4):
    #         cost += abs(identity_goal[i][j] - I4[i][j])

    return cost


def SU4_circuit(t_1, t_2, t_3, alpha_1, theta_1, beta_1, alpha_2, theta_2, beta_2, alpha_3, theta_3, beta_3, alpha_4, theta_4, beta_4):
    return np.kron(I2, R_z(beta_4)) @ np.kron(I2, R_y(theta_4)) @ np.kron(I2, R_z(alpha_4)) \
           @ np.kron(R_z(beta_3), I2) @ np.kron(R_y(theta_3), I2) @ np.kron(R_z(alpha_3), I2) \
           @ CNOT2 @ np.kron(I2, R_y(t_3)) @ CNOT1 @ np.kron(I2, R_y(t_2)) @ np.kron(R_z(t_1), I2) @ CNOT2 \
           @ np.kron(I2, R_z(beta_2)) @ np.kron(I2, R_y(theta_2)) @ np.kron(I2, R_z(alpha_2)) \
           @ np.kron(I2, R_z(beta_1)) @ np.kron(I2, R_y(theta_1)) @ np.kron(I2, R_z(alpha_1))


def cost_function_SU4(params: list):
    cost = 0
    SU4 = SU4_circuit(params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8],
                      params[9], params[10], params[11], params[12], params[13], params[14])
    identity_goal = SU4 @ np.linalg.inv(U)
    for i in range(4):
        for j in range(4):
            cost += abs(identity_goal[i][j] - I4[i][j])

    return cost

U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ CNOT1 @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2) @ CNOT2
# U = np.kron(pauli_X, I2) @ np.kron(I2, H) @ np.kron(I2, H) @ np.kron(I2, pauli_X) @ np.kron(pauli_X, I2)
# U = CNOT1 @ CNOT2 @ CNOT1

print(U)
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
        # parameters = np.array([np.pi / 7] * 6)
        parameters = np.array([np.pi / 7] * 15)
        # optimized_result = minimize(cost_function_SO4, parameters, method="Nelder-Mead")

        optimized_result = minimize(cost_function_SU4, parameters, method="L-BFGS-B")

        print("OPTIMIZED RESULT")
        print(optimized_result)
        print()

        # circuit_result = SO4_circuit(optimized_result['x'][0], optimized_result['x'][1], optimized_result['x'][2],
        #                              optimized_result['x'][3], optimized_result['x'][4], optimized_result['x'][5])
        circuit_result = SU4_circuit(optimized_result['x'][0], optimized_result['x'][1], optimized_result['x'][2],
                                     optimized_result['x'][3], optimized_result['x'][4], optimized_result['x'][5],
                                     optimized_result['x'][6], optimized_result['x'][7], optimized_result['x'][8],
                                     optimized_result['x'][9], optimized_result['x'][10], optimized_result['x'][11],
                                     optimized_result['x'][12], optimized_result['x'][13], optimized_result['x'][14])

rounded_circuit_result = round_circuit(circuit_result)
print("ROUNDED CIRCUIT RESULT")
print(rounded_circuit_result)

# print(circuit_result)
print()
print("U")
print(U)

print()

product = circuit_result @ np.linalg.inv(U)
rounded_product = round_circuit(product)
print("ROUNDED PRODUCT")
print(np.around(product, 2))
