import sympy as sp  # for the spar_to_dense helper function
import itertools as it
from flint import fmpq

class SparAlg:
    def __init__(self, ctx, d):
        self.d = d
        self.ctx = ctx
        self.zero = ctx.from_dict({})

    @staticmethod
    def cmul(x, y):
        a, b = x
        c, e = y
        return (a*c - b*e, a*e + b*c)

    @staticmethod
    def cadd_into(D, key, val, sign=1):
        if type(D) != dict:
            raise TypeError(f'D in cadd_into is not dict, but {type(D)}')
        if key in D.keys():
            re, im = D[key]    
            D[key] = (re + val[0], im + val[1]) if sign > 0 else (re - val[0], im - val[1])
        else:
            D[key] = val if sign > 0 else (-val[0], -val[1])
        return
    
    @staticmethod
    def clean(D):
        return {key: val for key, val in D.items() if not (val[0].is_zero() and val[1].is_zero())}

    def pad(self, A, left, right):
        if left == 0 and right == 0:
            return A
        d = self.d
        out = dict()
        for (r, c), v in A.items():
            for L in it.product(range(d), repeat=left):
                for R in it.product(range(d), repeat=right):
                    out[(L + r + R, L + c + R)] = v
        return out
    
    def spar_to_dense(self, A):
        d = self.d
        n = len(next(iter(A))[0])
        out = sp.zeros(d**n, d**n)
        for k, v in A.items():
            rtup, ctup = k
            n = len(rtup)
            r = sum([rtup[i] * d**(n-(i+1)) for i in range(n)])
            c = sum([ctup[i] * d**(n-(i+1)) for i in range(n)])
            out[r, c] = v
        return out

    # a bit suboptimal
    def mul(self, A, B):
        AB = dict()
        for ka, va in A.items():
            for kb, vb in B.items():
                if ka[1] == kb[0]:
                    prod = self.cmul(va, vb)
                    self.cadd_into(AB, (ka[0], kb[1]), prod)
        return self.clean(AB)

    def comm(self, A, B):
        out = self.mul(A, B)
        BA = self.mul(B, A)
        for k, v in BA.items():
            self.cadd_into(out, k, v, sign=-1)
        return self.clean(out)

    def add(self, A, B, sign=1):
        out = dict(A)
        for k, v in B.items():
            self.cadd_into(out, k, v, sign=sign)
        return out

    def scale(self, v, num, den):
        q = fmpq(num, den)
        return (v[0] * q, v[1]*q)

    # a, b could be inputs as well, if this solution turns out to be ugly
    def ti_comm(self, A, B, N):
        a = len(next(iter(A))[0])
        b = len(next(iter(B))[0])
        out = dict()
        for s in range(-(b-1), a):
            low = min(0, s)
            high = max(a, s+b)
            n = high - low
            Ap = self.pad(A, -low, high-a)
            Bp = self.pad(B, s-low, high-s-b)
            C = self.comm(Ap, Bp)
            Cp = self.pad(C, 0, N-n)
            for k, v in Cp.items():
                self.cadd_into(out, k, v, sign=1)
        return self.clean(out)

    def canonical(self, D, n):
        d = self.d
        cl_dict = dict()
        X_cur = D
        for l in range(n, 0, -1):
            if not X_cur:
                break
            r = dict()
            for (rw, cl), v in X_cur.items():
                # this is the sparse partial trace
                if rw[0] == cl[0]:
                    self.cadd_into(r, (rw[1:], cl[1:]), v)
            Y = dict(X_cur)
            for (rw, cl), v in r.items():
                for c in range(d):
                    self.cadd_into(Y, ((c,) + rw, (c,) + cl), self.scale(v, 1, d), -1)
            s = dict()
            for (rw, cl), v in Y.items():
                if rw[-1] == cl[-1]:
                    self.cadd_into(s, (rw[:-1], cl[:-1]), v)
            cl_op = dict(Y)
            for (rw, cl), v in s.items():
                for c in range(d):
                    self.cadd_into(cl_op, (rw + (c,), cl + (c,)), self.scale(v, 1, d), -1)
            cl_dict[l] = self.clean(cl_op)
            X_next = {}
            for (rw, cl), v in r.items():
                self.cadd_into(X_next, (rw, cl), self.scale(v, 1, d))
            for (rw, cl), v in s.items():
                self.cadd_into(X_next, (rw, cl), self.scale(v, 1, d))
            X_cur = self.clean(X_next) if l > 1 else {}
        return cl_dict