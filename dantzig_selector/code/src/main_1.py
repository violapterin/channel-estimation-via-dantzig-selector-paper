import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

import constants as cst
import classes as cls
import functions as fct

list_ss_nn_rr = [0.01 * 1.5 ** x for x in range (0,5)]
kk =fct.kk()
#for ss_nn_rr in list_ss_nn_rr:

pp = np.kron (
    fct.ff_bb().T @ fct.ff_rr().T @ kk.conj(),
    fct.ww_bb() @ fct.ww_rr() @ kk)
hh =fct.hh()
gg = kk.conj().T @ hh @ kk
g = fct.vectorize (gg)
z = fct.vectorize (fct.zz())

#y = pp @ g + z
y = pp @ g #XXX

gG = np.sqrt (2*np.log (cst.nn_h))
dd_ss = cls.Dd_ss (pp, y, gG);
dd_ss.run()

hh_hat = kk @ fct.inv_vectorize (dd_ss.g_hat, cst.nn_hh, cst.nn_hh) @ kk.conj().T
dd =hh_hat -hh

print(np.linalg.norm (hh, ord='fro'))
print(np.linalg.norm (hh_hat, ord='fro'))
print(np.linalg.norm (dd, ord='fro'))

#print(np.linalg.norm (ff_bb.val, ord=1))
#print(np.linalg.norm (ff_rr.val, ord=1))
#print(np.linalg.norm (ww_bb.val, ord=1))
#print(np.linalg.norm (ww_rr.val, ord=1))
#print(np.linalg.norm (hh.val, ord=1))
#print(np.linalg.norm (zz.val, ord=1))
#quit()
