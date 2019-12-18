#
# import numpy as np
# from qiskit.aqua import Operator
#
# from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, BasicAer
#
# q = QuantumRegister(1, "q")
# c = ClassicalRegister(1, "c")
# circ = QuantumCircuit(q, c)
#
# array = np.array([[0,1], [1,0]])
# cnot = Operator(matrix=array)
# evolve = cnot.evolve(evo_time=0, evo_mode="circuit", num_time_slices=1, quantum_registers=q)
# # print(evolve)
#
# circ += evolve
# circ.measure(q,c)
# print(circ)
#
# backend = BasicAer.get_backend("qasm_simulator")
#
# job = execute(circ, backend, shots=10)
# result = job.result()
# counts = result.get_counts(circ)
# print(counts)

# circ.h(q[0])
#
# circ.barrier()
# circ.measure(q, c)
#
# circ.draw(filename="circ.txt")
#
# backend = BasicAer.get_backend("qasm_simulator")
#
# job = execute(circ, backend, shots=10)
# result = job.result()
# counts = result.get_counts(circ)
#
# print(counts)



## personal testing

# theoretical_probs =  {'00000': 0.10139983459411607, '00001': 0.01735805702132643, '00010': 0.017358057021326437, '00011': 0.043801529890753046, '00100': 0.017358057021326423, '00101': 0.043801529890753046, '00110': 0.043801529890753046, '00111': 0.0017330570213264689, '01000': 0.017358057021326437, '01001': 0.04380152989075303, '01010': 0.043801529890753046, '01011': 0.0017330570213264678, '01100': 0.04380152989075303, '01101': 0.0017330570213264689, '01110': 0.0017330570213264678, '01111': 0.05942652989075301, '10000': 0.017358057021326437, '10001': 0.04380152989075307, '10010': 0.04380152989075308, '10011': 0.0017330570213264678, '10100': 0.04380152989075307, '10101': 0.0017330570213264689, '10110': 0.0017330570213264678, '10111': 0.05942652989075301, '11000': 0.04380152989075308, '11001': 0.0017330570213264695, '11010': 0.0017330570213264695, '11011': 0.05942652989075301, '11100': 0.0017330570213264654, '11101': 0.05942652989075301, '11110': 0.05942652989075301, '11111': 0.05933136172468943}
# experimental_probs = {'00000': 0.0997, '00001': 0.0175, '00010': 0.015, '00011': 0.0418, '00100': 0.017, '00101': 0.0459, '00110': 0.0443, '00111': 0.0018, '01000': 0.0155, '01001': 0.0471, '01010': 0.0451, '01011': 0.0015, '01100': 0.039, '01101': 0.0018, '01110': 0.0014, '01111': 0.0611, '10000': 0.0201, '10001': 0.0449, '10010': 0.043, '10011': 0.0021, '10100': 0.0446, '10101': 0.0023, '10110': 0.0022, '10111': 0.0582, '11000': 0.0435, '11001': 0.0018, '11010': 0.0016, '11011': 0.0595, '11100': 0.0024, '11101': 0.0631, '11110': 0.0612, '11111': 0.054}
#
#
# for key in theoretical_probs.keys():
#     if key in experimental_probs.keys():
#         print("present")
#     else:
#         print("not present")

from pyquil.quil import Gate
from pyquil.gates import H, CNOT, S
from pyquil import Program, get_qc
from pyquil.latex import to_latex
from pyquil.api import local_qvm

qvm = get_qc('3q-qvm')
program = Program(H(2))
program += S(2)
program += Gate.dagger(S(2))

print(to_latex(program))