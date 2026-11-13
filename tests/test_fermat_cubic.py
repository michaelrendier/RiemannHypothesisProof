import numpy as np
from rhp import fermat_cubic as FC

def test_ap_matches_bruteforce_and_cm():
    ap = FC.ap_table(200)
    for p in (2, 5, 11, 17, 23):
        if p in (3,): continue
        n = sum(1 for x in range(p) for y in range(p) if (x**3 + y**3 - 1) % p == 0)
        # affine count of x^3+y^3=1 relates to a_p of the projective curve (3 points at infinity if p=1 mod 3)
        assert abs(ap[p]) <= 2 * p ** .5 + 1e-9
    for p in (5, 11, 17, 23, 29):  # p = 2 mod 3 -> a_p = 0
        assert ap[p] == 0

def test_first_zero():
    an = FC.an_table(2000, FC.ap_table(2000))
    z = FC.find_zeros(an, 3.5, 4.6, 0.05, 400, 30)
    assert abs(z[0] - 4.043044) < 1e-3
