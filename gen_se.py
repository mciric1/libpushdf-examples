#!/usr/bin/env python3
#
#   Routines for generating worldsheet associahedron rational prefactors and the polynomial scattering equations
#

from itertools import combinations

from libpushdf import Monomial

def worldsheet_associahedron(n, Z):
    form_dnm = Z[0]
    prev = Z[0]

    for i in range(1,n-3):
        form_dnm *= (prev - Z[i])
        prev = Z[i]

    form_dnm *= (prev - 1)

    return -1 / form_dnm

def scattering_equations(n, RZ, RA):
    A_vars = RA.variables()

    def calc_a_var(i,j):
        # i,j are modulo, then sort so that i < j
        i = ((i - 1) % n) + 1
        j = ((j - 1) % n) + 1
        if i == j:
            return 0
        if i > j:
            i,j = j,i

        ind = (i - 1) * n + j - ((i + 2) * (i + 1) // 2) + (1 if i == 1 else 0)

        if j == i + 1 or (i,j) == (1,n):
            return 0

        assert ind <= len(A_vars)
        return A_vars[ind - 1]

    def calc_sA(A):
        sA = RA(0)
        for i,j in combinations(A,2):
            assert i < j
            sA += calc_a_var(i,j+1) + calc_a_var(i+1,j) - calc_a_var(i,j) - calc_a_var(i+1,j+1)
        return RZ(sA)

    def calc_zA(A):
        degrees = [0]*(n-3) # polynomial ring does not include z1, z[n-1] or z[n]
        for a in A:
            if a >= n - 1:
                continue # z[n-1] -> 1, z[n] -> infty
            degrees[a-1-1] = 1 # extra -1 for polynomial ring starting from z2
        zA = Monomial(RZ, degrees)
        return RZ(zA)

    # range starts from 2 since z[1] -> 0
    n_set = list(range(2, n + 1))

    P = []
    for l in range(1,(n-3)+1):
        P.append(sum(calc_zA(A) * calc_sA(A) for A in combinations(n_set,l+1) if n in A)) # condition for z[n] -> infty

    return P

# TODO migrate to unittest
if __name__ == '__main__':
    from libpushdf import QQ, PolynomialRing

    N = 6
    n_a_vars = N * (N - 3) // 2
    assert n_a_vars > 0

    RA = PolynomialRing(QQ, [f"a{i+1}" for i in range(n_a_vars)])
    RZ = PolynomialRing(RA.coeff_ring, [f"z{i+2}" for i in range(N-3)])
    RZ_A = PolynomialRing(RA, [f"z{i+2}" for i in range(N-3)])

    P = scattering_equations(N, RZ_A, RA)
    for p in P:
        print(p)

    print("")

    W = worldsheet_associahedron(N, RZ.variables())
    print(W)
