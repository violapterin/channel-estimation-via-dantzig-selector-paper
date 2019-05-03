import numpy as np
import cvxpy as cp
import constants as cst

class Dd_ss:
    def __init__(self, pp, y, gG):
        self.pp = pp
        self.y = y
        self.gG = gG
        self.g_hat = np.zeros ((cst.nn_hh))

    def run (self):
        for i in range (cst.nn_h):
            A.append (np.concatenate ([indication_repr_mat (i),
                                       np.zeros ((cst.nn_h, cst.nn_h))],
                                      axis= 1))
            b.append (np.zeros ((cst.nn_h)))
            c.append (np.concatenate ([np.zeros ((2 * cst.nn_h)), indication_vec (i)]))
            d.append (0)
        for i in range (cst.nn_h, 2 * cst.nn_h):
            A.append (np.concatenate ([-indication_repr_mat (i -cst.nn_h)
                                        @ find_repr_mat (self.pp.H)
                                        @ find_repr_mat (self.pp),
                                        np.zeros ((cst.nn_h, cst.nn_h))],
                                       axis= 1))
            b.append (indication_repr_mat (i -cst.nn_h)
                      @ find_repr_mat (self.pp.H)
                      @ find_repr_vec (self.y))
            c.append (np.zeros ((3 * cst.nn_h)))
            d.append (self.gG)
        t = np.concatenate ([np.zeros ((2 * cst.nn_h)), np.ones ((cst.nn_h))])
        x = cp.Variable (3 * cst.nn_h)

        constr = [cp.SOC (c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range (2 * cst.nn_h)]
        prob = cp.Problem (cp.Minimize (f.T@x), constr)
        prob.solve ()

        x_hat = x.value
        g_repr_hat = x_hat [0 : 2 * cst.nn_h -1]
        self.g_hat = inv_find_repr_vec (g_repr_hat)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Oo_mm_pp:
    def __init__(self):
        

class Ll_ss:
    def __init__(self, pp, y):




