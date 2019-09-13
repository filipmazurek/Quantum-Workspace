
import numpy as np
from qiskit.aqua import Operator

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, BasicAer

q = QuantumRegister(1, "q")
c = ClassicalRegister(1, "c")
circ = QuantumCircuit(q, c)

array = np.array([[0,1], [1,0]])
cnot = Operator(matrix=array)
evolve = cnot.evolve(evo_time=0, evo_mode="circuit", num_time_slices=1, quantum_registers=q)
# print(evolve)

circ += evolve
circ.measure(q,c)
print(circ)

backend = BasicAer.get_backend("qasm_simulator")

job = execute(circ, backend, shots=10)
result = job.result()
counts = result.get_counts(circ)
print(counts)

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



