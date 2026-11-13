import numpy as np
from rhp import chebyshev_bias as CB

def test_first_zeros_and_bias():
    g = CB.zeros(14.0)
    assert abs(g[0] - 6.0209) < 1e-3 and abs(g[1] - 10.2438) < 1e-3 and abs(g[2] - 12.9881) < 1e-3
    D = CB.D_exact(2000)
    assert abs(D[5] - (np.log(5) - np.log(3))) < 1e-12
