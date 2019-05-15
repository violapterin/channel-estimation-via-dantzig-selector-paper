import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

import constants as cst
import classes as cls
import functions as fct

list_ss_nn_rr = range()

for ss_nn_rr in list_ss_nn_rr:
    hh = cls.Hh()
    ff_bb = cls.Ff_bb()
    ff_rr = cls.Ff_rr()
    ww_bb = cls.Ww_bb()
    ww_rr = cls.Ww_rr()
    pp = np.kron ()


    dd_ss = cls.Dd_ss (pp, y, gG);
    hh_hat = cst.kk @ inv_vectorize (dd_ss.g_hat) @ cst.kk.conj().T

    oo_mm_pp = cls.Oo_mm_pp (pp, y);


