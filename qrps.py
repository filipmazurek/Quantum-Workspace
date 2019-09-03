"""
Filip Mazurek - 8/16/2019

An adaptation of Cat/Box/Scissors by Dr. James Wootton
https://medium.com/@decodoku/introducing-the-worlds-first-game-for-a-quantum-computer-50640e3c22e4

Made to be compatible with Qiskit (connecting to real quantum computer and local sim),
cirq simulation, pyquil simulation, and projectq simulation.

To be used as a simple reference for how to use each SDK.
"""

from qiskit.aqua.translators.ising import tsp
from qiskit.aqua.algorithms import VQE, ExactEigensolver

# print("\n\n\n\n===== Welcome to Cat/Box/Scissors! =====\n\n")
# print("  ~~ A game by the Decodoku project ~~ \n\n")
# print("When in doubt, press Enter key to continue!")
# input()
# print("You and your opponent choose one of two possible moves.")
# input()
# print("You win if your moves are the same.")
# input()
# print("Your opponent wins if they are different.")
# input()

print("How should this quantum program run?")

chosen = 0
runMethod = ""
while chosen == 0:
    runMethod = input("\nChoose your method (ibm_sim or ibm_real or cirq or pyquil or projq)\n")
    if (runMethod == "ibm_sim") | (runMethod == "ibm_real") | (runMethod == "cirq") \
            | (runMethod == "pyquil") | (runMethod == "pq") | (runMethod == "projq"):
        chosen = 1
    else:
        print("u wot m8? Try that again.")

# get human player to choose opponent
chosen = 0
opponent = -1
while chosen == 0:
    try:
        opponent = int(input("\nWhich qubit will be your opponent? (1,2,3, or 4)\n"))
        if (opponent >= 1) & (opponent <= 4):
            chosen = 1
        else:
            print("u wot m8? Try that again.")
    except ValueError:
        input("Try again")

# here 1 and 2 mean qubits 0 and 1, so change accordingly
# referee is always qubit 2
if opponent < 3:
    opponent = opponent - 1

# get human player to choose move
chosen = 0
humanMove = ""
while chosen == 0:
    humanMove = input("\nChoose your move (s or sdg)\n")
    if (humanMove == "s") | (humanMove == "sdg"):
        chosen = 1
    else:
        print("u wot m8? Try that again.")

# print("\nWe'll now send your move to the quantum referee - simulated.")
# input()
# print("It will take your opponents move and compare them.")
# input()

if (runMethod == "ibm_sim") | (runMethod == "ibm_real"):
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    from qiskit import execute
    from qiskit import BasicAer, IBMQ
    from qiskit.providers.ibmq import least_busy

    q = QuantumRegister(5)
    c = ClassicalRegister(1)
    qc = QuantumCircuit(q, c)

    if runMethod == "ibm_sim":
        backend = BasicAer.get_backend("qasm_simulator")
    else:
        import secrets
        provider = IBMQ.enable_account(secrets.IBM_TOKEN)
        large_enough_devices = provider.backends(filters=lambda x: x.configuration().n_qubits > 4 and
                                                               not x.configuration().simulator)
        backend = least_busy(large_enough_devices)
        print("The best IBM simulation backend is " + backend.name())

    # referee decides things in the X basis, so we prepare it in the |+> state
    qc.h(q[2])

    # implement human move
    if humanMove == "s":
        qc.s(q[2])
    else:
        qc.sdg(q[2])

    # opponent qubit is prepared in state |+> to randomly decide the move to make
    qc.h(q[opponent])

    # to implement the quantum move, first do an S in all cases
    qc.s(q[2])
    # then use a controlled-Z to make it into a Sdg if that's what the quantum player chooses
    qc.h(q[2])
    qc.cx(q[opponent], q[2])
    qc.h(q[2])

    # quantum player wins if the moves where different, which would leave referee in the |+> state
    # human player wins if the moves were the same
    # we measure in the X basis to see
    qc.h(q[2])
    qc.barrier()
    qc.measure(q[2], c)

    qc.draw(filename="circ.txt")

    # Run the quantum circuit on a statevector simulator backend
    job = execute(qc, backend, shots=1)
    result = job.result()
    counts = result.get_counts(qc)
    if "1" in counts:
        print("You win!")
    elif "0" in counts:
        print("You lose")
    else:
        print("Something weird happened")


if runMethod == "cirq":
    import cirq
    from cirq.ops import H, S, CNOT, measure

    qubits = [cirq.GridQubit(0, i) for i in range(5)]
    q2 = qubits[2]
    circuit = cirq.Circuit()

    circuit.append(H(q2))

    if humanMove == "s":
        circuit.append(S(q2))
    else:
        circuit.append(cirq.inverse(S(q2)))

    circuit.append([H(qubits[opponent]),
                    S(q2),
                    H(q2),
                    CNOT(qubits[opponent], q2),
                    H(q2),
                    H(q2),
                    measure(q2, key="m")
                    ])

    print(circuit)

    simulator = cirq.Simulator()
    results = simulator.run(circuit, repetitions=1)
    # print(results.histogram(key="m"))
    keys = results.histogram(key="m").keys()

    for key in keys:
        if key == 1:
            print("You win!")
        elif key == 0:
            print("You lose")
        else:
            print("Something weird happened")


# first, must start running compiler and vm using the following commands
# quilc -S
# qvm -S
if runMethod == "pyquil":
    from pyquil.quil import Gate
    from pyquil.gates import H, CNOT, S
    from pyquil import Program, get_qc
    from pyquil.api import local_qvm

    qvm = get_qc('5q-qvm')
    program = Program(H(2))

    if humanMove == "s":
        program += S(2)
    else:
        program += Gate.dagger(S(2))

    program += [H(opponent), S(2), H(2), CNOT(opponent, 2), H(2), H(2)]

    with local_qvm():
        results = qvm.run_and_measure(program, trials=1)

    if results[2][0] == 1:
        print("You win!")
    elif results[2][0] == 0:
        print("You lose")
    else:
        print("Something went weird")

    # print(results)

if runMethod == "projq":
    from projectq import MainEngine
    from projectq.ops import H, S, Sdag, CNOT, Measure

    eng = MainEngine()

    qubits = eng.allocate_qureg(5)
    q2 = qubits[2]

    H | q2

    if humanMove == "s":
        S | q2
    else:
        Sdag | q2

    H | qubits[opponent]
    S | q2
    H | q2
    CNOT | (qubits[opponent], q2)
    H | q2
    H | q2

    Measure | q2

    eng.flush()

    if int(q2) == 1:
        print("You win!")
    elif int(q2) == 0:
        print("You lose")
    else:
        print("Something went weird")