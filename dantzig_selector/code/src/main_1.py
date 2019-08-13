import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

import constants as cst
import classes as cls
import functions as fct

list_ss_nn_rr = [0.01 * 1.5 ** x for x in range (0,5)]
kk = fct.kk()

#for ss_nn_rr in list_ss_nn_rr:
ff_bb = cls.Ff_bb()
ff_rr = cls.Ff_rr()
ww_bb = cls.Ww_bb()
ww_rr = cls.Ww_rr()
hh = cls.Hh()
zz = cls.Zz()

ff_bb.gen()
ff_rr.gen()
ww_bb.gen()
ww_rr.gen()
hh.gen()
zz.gen()

pp = np.kron (ff_bb.val.T @ ff_rr.val.T @ kk.conj(), ww_bb.val @ ww_rr.val @ kk)
gg = kk.conj().T @ hh.val @ kk
g = fct.vectorize(gg)
z = fct.vectorize(zz.val)

y = pp @ g + 2*z

gG = np.sqrt (2*np.log (cst.nn_h))


dd_ss = cls.Dd_ss (pp, y, gG);
dd_ss.run()

hh_hat = kk @ fct.inv_vectorize (dd_ss.g_hat, cst.nn_hh, cst.nn_hh) @ kk.conj().T

print(np.linalg.norm (hh_hat, ord=1))

#print(np.linalg.norm (ff_bb.val, ord=1))
#print(np.linalg.norm (ff_rr.val, ord=1))
#print(np.linalg.norm (ww_bb.val, ord=1))
#print(np.linalg.norm (ww_rr.val, ord=1))
#print(np.linalg.norm (hh.val, ord=1))
#print(np.linalg.norm (zz.val, ord=1))
#quit()
