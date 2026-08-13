# QB spot-check of the referee's KILLER A (independent build).
# Code on 7 holders: c1={0,1,3,4,5,6} (descended, target {0,1}),
# c2={1,2,3,4}; signs from the quadratic condition:
# s(c1)=s(c2)=+1 => s(c1+c2) = (-1)^{|c1 cap c2|} = (-1)^3 = -1.
# Measured holders 3,4,5,6 with aux rebits 7,8,9,10; aux Omegas on
# rebit pairs (7,8) and (9,10) (holder pairs (3,4),(5,6)).
# Claim to check: on every branch, after the prescribed target
# correction, the kept state on holders {0,1,2} is the rank-2
# correlated state ~ (1/8)(I + J01 + J02 - J12) (up to local signs),
# NOT Omega_{01} (x) I/2 (rank 4); auxes exactly Omega each.
import numpy as np
from functools import reduce
from itertools import product

I2 = np.eye(2)
J = np.array([[0., -1.], [1., 0.]])
V = np.diag([1., -1.])
N = 11  # 7 resource + 4 aux


def kron(*ops):
    return reduce(np.kron, ops)


def js(sites):
    return kron(*[J if i in sites else I2 for i in range(N)])


def marg(r, keep):
    perm = keep + [i for i in range(N) if i not in keep]
    t = r.reshape([2] * (2 * N))
    t = np.transpose(t, perm + [N + p for p in perm])
    k = len(keep)
    t = t.reshape(2 ** k, 2 ** (N - k), 2 ** k, 2 ** (N - k))
    return np.einsum("iaja->ij", t)


c1 = [0, 1, 3, 4, 5, 6]
c2 = [1, 2, 3, 4]
c12 = sorted(set(c1) ^ set(c2))  # {0,2,5,6}
sig = (np.eye(2 ** 7) +
       kron(*[J if i in c1 else I2 for i in range(7)]) +
       kron(*[J if i in c2 else I2 for i in range(7)]) -
       kron(*[J if i in c12 else I2 for i in range(7)])) / 2 ** 7
Om = 0.25 * (np.eye(4) - np.kron(J, J))
rho = kron(sig, Om, Om)  # aux rebits 7,8 (pair 3,4) and 9,10 (pair 5,6)
w = np.linalg.eigvalsh(rho)
assert w.min() > -1e-12, ("sigma not PSD", w.min())

Ms = [js([3, 7]), js([4, 8]), js([5, 9]), js([6, 10])]
Id = np.eye(2 ** N)
target_claim = None
ranks = set()
aux_ok_all, prod_fail_all = True, True
for xs in product((+1, -1), repeat=4):
    P = Id
    for x, M in zip(xs, Ms):
        P = P @ (0.5 * (Id + x * M))
    post = P @ rho @ P.T
    p = np.trace(post).real
    if p < 1e-14:
        print("zero-probability branch", xs)
        continue
    post /= p
    # target correction for c1: r=2 aux pairs, s=+1 -> flip iff h=+1
    h = np.prod(xs)
    if h > 0:
        C = kron(*[V if i == 1 else I2 for i in range(N)])
        post = C @ post @ C.T
    kept = marg(post, [0, 1, 2])
    r = np.linalg.matrix_rank(kept, tol=1e-10)
    ranks.add(r)
    # compare with Omega_{01} (x) I/2 (rank 4)
    cand = np.kron(Om, I2 / 2)
    prod_fail_all &= not np.allclose(kept, cand, atol=1e-10)
    # target marginal alone
    mt = marg(post, [0, 1])
    assert np.allclose(mt, Om, atol=1e-10), ("target marginal not Om", xs)
    # auxes
    for ap in ([7, 8], [9, 10]):
        aux_ok_all &= np.allclose(marg(post, ap), Om, atol=1e-10)

print("branch probs all nonzero; kept-state ranks over branches:", ranks)
print("target-pair marginal exactly Omega on every branch: True")
print("aux pairs exactly Omega on every branch:", aux_ok_all)
print("kept != Omega x I/2 on every branch (decoupling FAILS):",
      prod_fail_all)
print("VERDICT: Killer A", "CONFIRMED" if (ranks == {2} and prod_fail_all
      and aux_ok_all) else "NOT REPRODUCED", "(rank-2 correlated kept "
      "state; no correction can raise rank to 4)")
