from functools import reduce
import itertools
import numpy as np
from qiskit.quantum_info import Pauli
from qiskit.aqua import Operator

p = 10
beta = np.random.uniform(0, np.pi*2, p)
gamma = np.random.uniform(0, np.pi*2, p)

print(beta)