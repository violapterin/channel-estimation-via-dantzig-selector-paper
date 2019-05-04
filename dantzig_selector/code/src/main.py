import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

import constants as cst
import classes as cls
import functions as fct

list_ss_nn_rr = range()

for ss_nn_rr in list_ss_nn_rr:
    hh =Hh()
    yy =hh
    dd_ss = Dd_ss (pp, y, gG);
    hh_hat = cst.kk @ inv_vectorize (dd_ss.g_hat) @ cst.kk.conj().T

    oo_mm_pp = Oo_mm_pp (pp, y);


