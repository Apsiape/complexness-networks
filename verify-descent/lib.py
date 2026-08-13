"""Frame-code descent referee.  Two INDEPENDENT engines.

Engine A : dense real matrices (numpy), no use of the commuting structure.
Engine B : symbolic J-string coefficient algebra.

Conventions
-----------
rebit = R^2 ;  J = [[0,-1],[1,0]] ;  J^2 = -I ; J^T = -J.
Sites 0..n-1.  Site k occupies bit (n-1-k) of the computational index
(standard kron(A_0,...,A_{n-1}) ordering).
J_S |x> = (-1)^{popcount(x & S)} |x xor S|.
"""
import numpy as np
import itertools

# ---------------------------------------------------------------- basics
Jl = np.array([[0.0, -1.0], [1.0, 0.0]])
Il = np.eye(2)
Zl = np.diag([1.0, -1.0])


def mask_of(sites, n):
    m = 0
    for k in sites:
        m |= 1 << (n - 1 - k)
    return m


def popcount(a):
    a = np.asarray(a, dtype=np.int64).copy()
    c = np.zeros_like(a)
    while a.any():
        c += a & 1
        a >>= 1
    return c


def jstring_dense(sites, n):
    """dense J_S built by index arithmetic (cheap, exact)."""
    N = 1 << n
    m = mask_of(sites, n)
    cols = np.arange(N, dtype=np.int64)
    rows = cols ^ m
    vals = np.where(popcount(cols & m) & 1, -1.0, 1.0)
    M = np.zeros((N, N))
    M[rows, cols] = vals
    return M


def jstring_kron(sites, n):
    """dense J_S built by literal kron -- used only to validate jstring_dense."""
    out = np.array([[1.0]])
    for k in range(n):
        out = np.kron(out, Jl if k in sites else Il)
    return out


def ptrace_dense(rho, n, traced):
    """trace out the listed sites; returns operator on the remaining sites
    in their original relative order."""
    cur = rho
    sites = list(range(n))
    for t in sorted(traced, reverse=True):
        k = sites.index(t)
        nn = len(sites)
        d1 = 1 << k
        d2 = 1 << (nn - k - 1)
        T = cur.reshape(d1, 2, d2, d1, 2, d2)
        cur = np.einsum('aibcid->abcd', T).reshape(d1 * d2, d1 * d2)
        sites.pop(k)
    return cur, sites


def local_Z(site, n):
    out = np.array([[1.0]])
    for k in range(n):
        out = np.kron(out, Zl if k == site else Il)
    return out


# ------------------------------------------------- Engine B : symbolic
class Sym(dict):
    """operator = dict frozenset(sites) -> coefficient."""
    pass


def sym_scale(A, c):
    return Sym({k: v * c for k, v in A.items()})


def sym_add(A, B):
    out = Sym(A)
    for k, v in B.items():
        out[k] = out.get(k, 0.0) + v
        if abs(out[k]) < 1e-14:
            del out[k]
    return out


def sym_mul(A, B):
    """J_S J_T = (-1)^{|S cap T|} J_{S xor T}"""
    out = Sym()
    for S, a in A.items():
        for T, b in B.items():
            sgn = -1.0 if (len(S & T) % 2) else 1.0
            K = S ^ T
            out[K] = out.get(K, 0.0) + sgn * a * b
    return Sym({k: v for k, v in out.items() if abs(v) > 1e-14})


def sym_ptrace(A, traced):
    """Tr_{traced}: strings touching a traced site vanish; others pick 2^{|traced|}."""
    f = float(1 << len(traced))
    tr = frozenset(traced)
    return Sym({S: v * f for S, v in A.items() if not (S & tr)})


def sym_trace(A, n):
    return A.get(frozenset(), 0.0) * float(1 << n)


def sym_to_dense(A, n, sites_order=None):
    """sites_order: list of global site labels in the order they occupy the
    tensor factors of the dense matrix."""
    if sites_order is None:
        sites_order = list(range(n))
    idx = {s: i for i, s in enumerate(sites_order)}
    nn = len(sites_order)
    M = np.zeros((1 << nn, 1 << nn))
    for S, v in A.items():
        M += v * jstring_dense([idx[s] for s in S], nn)
    return M


def sym_I():
    return Sym({frozenset(): 1.0})


def sym_J(sites):
    return Sym({frozenset(sites): 1.0})


# ------------------------------------------------------- code helpers
def code_words(gens):
    """gens: list of frozensets. returns list of (word_frozenset, subset_tuple)"""
    out = []
    k = len(gens)
    for bits in itertools.product([0, 1], repeat=k):
        w = frozenset()
        for i, b in enumerate(bits):
            if b:
                w = w ^ gens[i]
        out.append((w, bits))
    return out


def sign_from_q(bits, gens, qvals):
    """s_c = (-1)^{q(c)},  q(sum_{i in T} g_i) = sum q(g_i) + sum_{i<j in T} g_i . g_j
    guarantees s_u s_v (-1)^{u.v} = s_{u+v}."""
    T = [i for i, b in enumerate(bits) if b]
    q = 0
    for i in T:
        q += qvals[i]
    for a in range(len(T)):
        for b in range(a + 1, len(T)):
            q += len(gens[T[a]] & gens[T[b]])
    return -1.0 if (q % 2) else 1.0
