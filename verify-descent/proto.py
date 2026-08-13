"""Protocol runner: two independent engines + the closed-form prediction."""
import numpy as np, itertools, sys
from lib import *


class Setup:
    """holders 0..H-1 (one resource rebit each, site index = holder).
       aux pairs: list of (u,v) holder pairs sharing an Omega; aux rebits get
       site indices H + 2t (for u) and H + 2t + 1 (for v).
       measured holders U = all holders in aux pairs.
       traced_extra = holders whose resource rebit is also discarded (W)."""

    def __init__(self, H, auxpairs, gens, qvals=None, traced_extra=()):
        self.H = H
        self.auxpairs = [tuple(p) for p in auxpairs]
        self.P = len(self.auxpairs)
        self.n = H + 2 * self.P
        self.aux = {}                      # holder -> aux site
        for t, (u, v) in enumerate(self.auxpairs):
            self.aux[u] = H + 2 * t
            self.aux[v] = H + 2 * t + 1
        self.U = sorted(self.aux.keys())
        self.W = sorted(traced_extra)
        self.V = [h for h in range(H) if h not in self.U and h not in self.W]
        self.traced = sorted(self.U + self.W)
        self.keep = [s for s in range(self.n) if s not in self.traced]
        self.gens = [frozenset(g) for g in gens]
        self.qvals = list(qvals) if qvals is not None else [0] * len(self.gens)
        self.words = []                    # (word, bits, sign)
        for w, bits in code_words(self.gens):
            self.words.append((w, bits, sign_from_q(bits, self.gens, self.qvals)))
        # sanity: even code
        assert all(len(w) % 2 == 0 for w, _, _ in self.words), "code has odd-weight word"

    def auxpair_index(self, holder):
        for t, (u, v) in enumerate(self.auxpairs):
            if holder in (u, v): return t
        return None


# ------------------------------------------------------------ Engine A
def engineA(S, x):
    """dense.  x: dict holder->+-1 for holders in S.U .
       returns (prob, rho_out_dense, keep_sites)"""
    n = S.n
    N = 1 << n
    rho = np.zeros((N, N))
    for w, bits, s in S.words:
        rho += s * jstring_dense(sorted(w), n)
    rho *= 2.0 ** (-S.H)
    # tensor the Omegas (they live on aux sites; multiply in)
    for (u, v) in S.auxpairs:
        a, b = S.aux[u], S.aux[v]
        Om = 0.25 * (np.eye(N) - jstring_dense([a, b], n))
        rho = rho @ Om
    # projector
    Pi = np.eye(N)
    for u in S.U:
        M = jstring_dense([u, S.aux[u]], n)
        Pi = Pi @ (0.5 * (np.eye(N) + x[u] * M))
    cond = Pi @ rho @ Pi
    prob = float(np.trace(cond))
    out, keep = ptrace_dense(cond, n, S.traced)
    return prob, out, keep


# ------------------------------------------------------------ Engine B
def engineB(S, x):
    """symbolic J-string algebra."""
    rho = Sym()
    for w, bits, s in S.words:
        rho[frozenset(w)] = rho.get(frozenset(w), 0.0) + s * 2.0 ** (-S.H)
    for (u, v) in S.auxpairs:
        a, b = S.aux[u], S.aux[v]
        Om = sym_add(sym_scale(sym_I(), 0.25), sym_scale(sym_J([a, b]), -0.25))
        rho = sym_mul(rho, Om)
    Pi = sym_I()
    for u in S.U:
        Pu = sym_add(sym_scale(sym_I(), 0.5), sym_scale(sym_J([u, S.aux[u]]), 0.5 * x[u]))
        Pi = sym_mul(Pi, Pu)
    cond = sym_mul(sym_mul(Pi, rho), Pi)
    prob = sym_trace(cond, S.n)
    out = sym_ptrace(cond, S.traced)
    return prob, out


# ------------------------------------------------- closed-form prediction
def predict(S, x, verbose=False):
    """The claimed analytic form.
       For each codeword c':  D = c' cap U, R = c' \\ U.
       covered  iff  R subset V (i.e. c' misses W).
       aligned  iff  D is a union of complete aux pairs.
       r        = # aux pairs fully inside D.
       coeff    = s_c * (-1)^{|D|} * (-1)^{r} * h,   h = prod_{u in D} x_u.
       Returns (sym operator on kept sites, diagnostics)."""
    out = Sym()
    diag = []
    all_aligned = True
    for w, bits, s in S.words:
        D = frozenset(w) & frozenset(S.U)
        R = frozenset(w) - frozenset(S.U)
        covered = not (R & frozenset(S.W))
        aligned = True
        r = 0
        for t, (u, v) in enumerate(S.auxpairs):
            inD = (u in D) + (v in D)
            if inD == 2: r += 1
            elif inD == 1: aligned = False
        h = 1.0
        for u in D: h *= x[u]
        coeff = s * ((-1.0) ** len(D)) * ((-1.0) ** r) * h
        diag.append(dict(word=sorted(w), s=s, D=sorted(D), R=sorted(R), covered=covered,
                         aligned=aligned, r=r, h=h, coeff=coeff))
        if not aligned: all_aligned = False
        if covered and aligned:
            key = frozenset(R)
            out[key] = out.get(key, 0.0) + coeff
    # tensor the untouched Omegas
    res = Sym({k: v for k, v in out.items() if abs(v) > 1e-14})
    for (u, v) in S.auxpairs:
        a, b = S.aux[u], S.aux[v]
        Om = sym_add(sym_scale(sym_I(), 0.5), sym_scale(sym_J([a, b]), -0.5))
        res = sym_mul(res, Om)
    if verbose:
        for d in diag:
            print("   c'=%-22s s=%+d D=%-14s R=%-10s cov=%d align=%d r=%d h=%+d coeff=%+d"
                  % (d['word'], d['s'], d['D'], d['R'], d['covered'], d['aligned'],
                     d['r'], d['h'], d['coeff']))
    return res, diag, all_aligned


# ------------------------------------------------------------ utilities
def normalise(sym_op, nkeep):
    tr = sym_trace(sym_op, nkeep)
    if abs(tr) < 1e-13: return None
    return sym_scale(sym_op, 1.0 / tr)


def dense_to_sym(M, keep):
    """extract J-string coefficients of a dense operator on |keep| sites."""
    nn = len(keep)
    out = Sym()
    for r in range(nn + 1):
        for S in itertools.combinations(range(nn), r):
            c = np.trace(jstring_dense(S, nn) @ M) / (((-1.0) ** r) * (1 << nn))
            if abs(c) > 1e-10:
                out[frozenset(keep[i] for i in S)] = c
    return out


def support_group(sym_op):
    return set(k for k, v in sym_op.items() if abs(v) > 1e-10)


def is_disjoint_pair_generated(G):
    """G: set of frozensets forming a group under xor. True iff generated by
    pairwise-disjoint weight-2 elements."""
    pairs = [g for g in G if len(g) == 2]
    # pairwise disjoint?
    for a, b in itertools.combinations(pairs, 2):
        if a & b: return False, "pairs intersect: %s %s" % (sorted(a), sorted(b))
    # do they generate G?
    gen = {frozenset()}
    for p in pairs:
        gen = gen | {g ^ p for g in gen}
        # close
        changed = True
        while changed:
            changed = False
            new = {a ^ b for a in gen for b in gen}
            if not new <= gen:
                gen |= new; changed = True
    if gen != G:
        return False, "pairs do not generate the surviving group"
    return True, "ok"
