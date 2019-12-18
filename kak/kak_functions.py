import numpy as np
import math
from typing import Tuple
import cmath


abs_error = pow(10, -8)
magic = 1/math.sqrt(2) * np.array([[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]])
eps = np.finfo(float).eps
sigx = np.array([[0, 1], [1, 0]])
sigy = np.array([[0, -1j], [1j, 0]])
sigz = np.array([[1, 0], [0, -1]])
sigxx = np.kron(sigx, sigx)
sigyy = np.kron(sigy, sigy)
sigzz = np.kron(sigz, sigz)

tester = np.array([[1, 0, 0, 2], [0, 2, 1, 0], [0, 2, -1, 0], [1, 0, 0, -2]])


def kak1(U: np.array) -> Tuple[np.array, np.array, np.array, np.array, np.array]:
    # Calculates a special case (KAK1)
    # of Cartan's KAK Decomposition.
    # For an input 4d unitary matrix U,
    # it finds [A1, A0, class_vec, B1, B0]
    # such that U = kron(A1,A0)*exp(i*M)*kron(B1,B0)
    # where matrix M is a function of
    # class_vec(1:4)

    # If
    # k=class_vec
    # and
    # M=k(4)+sigxx*k(1)+sigyy*k(2)+sigzz*k(3)
    # then
    # kron(A1,A0)*exp(i*M)*kron(B1,B0) = U

    if np.shape(U)[0] != 4 or np.shape(U)[1] != 4:
        raise ValueError("input matrix for kak decomposition is not 4x4")
    if np.linalg.norm(U @ U.conjugate().transpose() - np.eye(4)) > 1e-8:
        raise ValueError("input matrix for kak decomposition is not unitary")

    # normalize U so that it has
    # unit determinant. Store
    # det(U)^(1/4) for later use

    root4th = pow(np.linalg.det(U), (1 / 4))

    if abs(root4th - 1) > 1e-8:
        U = U / root4th

    U_bell = magic.conjugate().transpose() @ U @ magic

    [L, D, R] = odo(U_bell)

    # print("L\n", L, "\nD\n", D, "\nR\n", R)

    class_vec = diag_to_class_vec(np.diagonal(D))
    class_vec[3] = class_vec[3] + np.angle(root4th)

    [A1, A0] = SO4_to_SU2xSU2(L)
    [B1, B0] = SO4_to_SU2xSU2(R.conjugate().transpose())

    return A1, A0, class_vec, B1, B0


def class_k_to_mat(k: np.array) -> np.array:
    # Builds a 4d matrix mat = e^{i (k_0 + k.Sigma)}
    # characterized by a class vector k(1:3) and
    # a global phase k(4).

    mat = np.exp(1j * k[3]) \
          @ (np.cos(k[0]) * np.eye(4) + 1j * np.sin(k[0]) * sigxx) \
          @ (np.cos(k[1]) * np.eye(4) + 1j * np.sin(k[1]) * sigyy) \
          @ (np.cos(k[2]) * np.eye(4) + 1j * np.sin(k[2]) * sigzz)

    return mat


def canonical_class_vec(a: np.array) -> np.array:
    # Reduces a class vector to its canonical form.
    #
    # Input "a" is 4d real column vector, a(1:3) is a class vector and  a(4) is a global phase.
    #
    # Output "b" is the new 4d real column vector,
    # b(1:3) is a canonical vector equivalent to a(1:3) and b(4)=a(4).
    # b will satisfy:
    # pi/2 > bx >= by >= bz >= 0
    # pi/2 >= bx + by
    # if bz==0, then pi/4 >= bx

    b = a
    pih = np.pi / 2

    for j in range(2):
        while b[j] < 0:
            b[j] = b[j] + pih
        while b[j] >= pih:
            b[j] = b[j] - pih

    b[0:2] = np.flipud(np.sort(b[0: 2]))  # now b in increasing order
    # a good example for the raison d'etre of the next if() is
    # b1=pi/2-eps, b2=pi/2-2eps, b3=pi/2-3eps,
    # where eps>0 and small

    if b[0] + b[1] > pih:
        old_b0 = b[0]
        b[0] = pih - b[1]
        b[1] = pih - old_b0
        # at this point b(1)>b(2) but b(3) may be larger than b(1) or b(2)
        if b[2] > b[1]:
            b[0:2] = np.flipud(np.sort(b[0:2]))

    if b[2] < 1e-14 and b[0] > np.pi / 4:
        b[0] = pih - b[0]

    return b


def diag_to_class_vec(dd: np.array) -> np.array:
    # Finds a class vector cvec(1:3) and global phase cvec(4) from the
    # diagonal entries dd(1:4) of a diagonal unitary matrix.
    t1 = np.angle(dd[0])
    t2 = np.angle(dd[1])
    t3 = np.angle(dd[2])
    t4 = np.angle(dd[3])

    k0 = (t1 + t2 + t3 + t4) / 4
    kx = (t1 + t2 - t3 - t4) / 4
    ky = (-t1 + t2 - t3 + t4) / 4
    kz = (t1 - t2 - t3 + t4) / 4

    cvec = np.array([kx, ky, kz, k0]).transpose()

    return cvec


def rjd(A: np.array, threshold=np.sqrt(eps)) -> Tuple[np.array, np.array]:
    # TODO: does this make sense?
    m = np.shape(A)[0]
    nm = np.shape(A)[1]
    A = np.array(A)
    V = np.eye(m)

    encore = True

    # TODO: check if this while loop matches
    while encore:
        encore = False
        for p in range(0, m-1):
            for q in range(p+1, m):
                top_g = np.array(A[p, p:nm:m] - A[q, q:nm:m])
                bot_g = np.array(A[p, q:nm:m] + A[q, p:nm:m])
                g = np.vstack([top_g, bot_g])
                # g = np.array([[A[p, p:nm:m] - A[q, q:nm:m]], [A[p, q:nm:m] + A[q, p:nm:m]]])
                g = g @ g.conjugate().transpose()
                ton = g[0, 0] - g[1, 1]
                toff = g[0, 1] + g[1, 0]
                theta = 0.5 * cmath.atan(toff / (ton + np.sqrt(ton * ton + toff * toff)))
                c = np.cos(theta)
                s = np.sin(theta)
                encore = encore or (abs(s) > threshold)

                if abs(s) > threshold:
                    Mp = A[:, p:nm:m]
                    Mq = A[:, q:nm:m]
                    A[:, p:nm:m] = c * Mp + s * Mq
                    A[:, q:nm:m] = c * Mq - s * Mp
                    rowp = A[p, :]
                    rowq = A[q, :]
                    A[p, :] = c * rowp + s * rowq
                    A[q, :] = c * rowq - s * rowp
                    temp = V[:, p]
                    V[:, p] = c * V[:, p] + s * V[:, q]
                    V[:, q] = c * V[:, q] - s * temp

    qDs = A

    return V, qDs


def simul_real_svd(X: np.array, Y: np.array) -> Tuple[np.array, np.array, np.array, np.array]:
    dim_U = np.shape(X)[0]

    if (np.shape(X)[0] != np.shape(X)[0]) or (np.shape(X)[1] != np.shape(Y)[1]):
        raise ValueError("X and Y don't have the same dimensions")
    print("\nX\n", X)
    [LX, DX, RX] = np.linalg.svd(X)
    DX = np.diag(DX)  # converting to diagonal matrix

    # print("\nLX\n", LX, "\nDX\n", DX, "\nRX\n", RX)

    tol = dim_U * DX[0, 0] * eps
    # rank_X = np.linalg.matrix_rank(X)
    # print("\nrank_X\n", rank_X)
    rank_X = 1
    j = 2 - 1

    done = False
    while j < dim_U and not done:
        if DX[j, j] > tol:
            rank_X += 1
        else:
            done = True
        j += 1

    # print("\nrank_X\n", rank_X)

    if rank_X < dim_U:
        DX00 = DX[:rank_X, :rank_X]
        Z01 = np.zeros(rank_X, dim_U - rank_X)
        Z10 = np.zeros(dim_U - rank_X, rank_X)

        G = LX.conjugate().transpose() @ Y @ RX
        G00 = G[:rank_X, :rank_X]
        G11 = G[rank_X:dim_U:1, rank_X:dim_U:1]

        [P00, qDs] = rjd([DX00, G00], 1e-10)
        DG00 = qDs[:rank_X - 1, rank_X:2*rank_X:1]

        [LG11, DG11, RG11] = np.linalg.svd(G11)
        DY = [[DG00, Z01], [Z10, DG11]]

        R = RX @ np.array([[P00, Z01], [Z10, RG11]])
        L = LX @ [[P00, Z01], [Z10, LG11]]

    else:
        DX00 = DX
        G00 = LX.conjugate().transpose() @ Y @ RX
        [P00, qDs] = rjd([DX00, G00], 1e-12)
        qDs = np.array(qDs)
        DY = qDs[:dim_U, dim_U:(2*dim_U):1]
        R = RX @ P00
        L = LX @ P00

    # TODO: potential bug. Can use isclose if checking for det = -1
    if np.linalg.det(L) < 0:
        L[:, 0] = -L[:, 0]
        DX[0, 0] = -DX[0, 0]
        DY[0, 0] = -DY[0, 0]

    if np.linalg.det(R) < 0:
        R[:, 0] = -R[0, 0]
        DY[0, 0] = -DY[0, 0]

    return L, DX, DY, R


def SO4_to_SU2xSU2(Q: np.array) -> Tuple[np.array, np.array]:
    # Given a matrix Q\in SO(4), this function
    # finds matrices A, B\in SU(2) such that Q = M' * kron(A, B) * M

    if np.linalg.norm(Q * Q.transpose() - np.eye(4)) > abs_error:
        raise ValueError("input matrix for SO4_to_SU2xSU2 is not orthogonal")

    if abs(np.linalg.det(Q) - 1) > abs_error:
        raise ValueError("input matrix for SO4_to_SU2xSU2 is not special")

    U = magic @ Q @ magic.conjugate().transpose()

    # Start by assuming that U = [ a11  * B, a12  * B]  Will check validity of this assumption at the end
    #                            [-a12' * B, a11' * B]

    A = np.zeros((2, 2))
    B = np.zeros((2, 2))

    a11_sq_mat = U[:2, :2] @ U[2:, 2:].conjugate().transpose()
    a11_sq = a11_sq_mat[0, 0]

    a12_sq_mat = -1 * (U[:2, 2:] @ U[2:, :2].conjugate().transpose())
    a12_sq = a12_sq_mat[0, 0]

    a11_a12H_mat = U[:2, :2] @ U[:2, 2:].conjugate().transpose()
    a11_a12H = a11_a12H_mat[0, 0]

    A[0, 0] = math.sqrt(a11_sq)
    A[0, 1] = math.sqrt(a12_sq)

    # sqrt ambiguous in its sign, so the values just assigned to A(1,1) and to A(1,2) may have the wrong sign.
    # Fix relative sign:

    if abs(A[0, 0] * A[0, 1].conjugate().transpose() - a11_a12H) > abs_error:
        A[0, 1] = -A[0, 1]

    A[1, 0] = -A[0, 1].conjugate().transpose()
    A[1, 1] = A[0, 0].conjugate().transpose()

    if abs(A[0, 0]) > abs(A[0, 1]):
        B = U[:2, :2] / A[0, 0]
    else:
        B = U[:2, 2:] / A[0, 1]

    if np.linalg.norm(np.kron(A, B) - U) > abs_error:
        raise ValueError("The error in  SO4_to_SU2xSU2 is too large")

    return A, B


def odo(U: np.array) -> Tuple[np.array, np.array, np.array, np.array]:
    # odo=orthogonal diagonal orthogonal
    # Given a  unitary matrix U, it finds special (det=1) orthogonal matrices L and R and a diagonal unitary matrix D
    # such that L * D * R' = U

    dim_U = np.shape(U)[0]
    if np.linalg.norm(U @ U.conjugate().transpose() - np.eye(dim_U)) > abs_error:
        raise ValueError("input matrix for odo decomposition is not unitary")

    X = (U + U.conjugate()) / 2
    Y = (U - U.conjugate()) / 2j
    # print("\nX\n", X, "\nY\n", Y)
    [L, DX, DY, R] = simul_real_svd(X, Y)

    return L, DX, DY, R


pauli_X = np.array([[0, 1], [1, 0]])
I2 = np.eye(2)
I4 = np.eye(4)
in_U = np.kron(pauli_X, I2)
# print(in_U)
in_U2 = np.array([[0,0,0,1],[-1,0,0,0],[0,1,0,0], [0,0,1,0]])
# print(in_U)
# test_mat = np.array([[0,1,2,3], [4,5,6,7],[8,9,10,11],[12,13,14,15]])
# print(test_mat)
# print()
# print(test_mat[0:4:2, 0:4:1])

[A1, A0, class_vec, B1, B0] = kak1(in_U2)
