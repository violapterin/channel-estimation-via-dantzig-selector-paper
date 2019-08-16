import numpy as np
import cvxpy as cp

import constants as cst
import functions as fct


class Dd_ss:
    def __init__(self, pp, y, gG):
        self.pp = pp
        self.y = y
        self.gG = gG
        self.g_hat = np.zeros ((cst.nn_hh))

    def run (self):
        aa=[]
        b=[]
        c=[]
        d=[]
        for i in range (cst.nn_h):
            aa.append (
                np.concatenate (
                    [fct.indication_repr_mat (i),
                        np.zeros ((2 * cst.nn_h, cst.nn_h))],
                    axis= 1))
            b.append (np.zeros ((2 * cst.nn_h)))
            c.append (
                np.concatenate (
                    [np.zeros ((2 * cst.nn_h)), fct.indication_vec (i)]))
            d.append (0)
        for i in range (cst.nn_h):
            aa.append (
                np.concatenate (
                    [-fct.indication_repr_mat (i)
                            @ fct.find_repr_mat (self.pp.conj().T)
                            @ fct.find_repr_mat (self.pp),
                        np.zeros ((2 * cst.nn_h, cst.nn_h))],
                    axis= 1))
            b.append (fct.indication_repr_mat (i)
                @ fct.find_repr_mat (self.pp.conj().T)
                @ fct.find_repr_vec (self.y))
            c.append (np.zeros ((3 * cst.nn_h)))
            d.append (self.gG)
        t = np.concatenate (
            [np.zeros ((2 * cst.nn_h)),
                np.ones ((cst.nn_h))])
        x = cp.Variable (3 * cst.nn_h)

        constr = [cp.SOC (c[i].T @ x + d[i], aa[i] @ x + b[i]) for i in range (2 * cst.nn_h)]
        prob = cp.Problem (cp.Minimize (t.T@x), constr)
        prob.solve ()

        x_hat = x.value
        g_repr_hat = x_hat [0 : 2 * cst.nn_h]
        self.g_hat = fct.inv_find_repr_vec (g_repr_hat)

# TODO
class Oo_mm_pp:
    def __init__(self, pp, y, eG):
        self.pp = pp
        self.y = y
        self.g_hat = np.zeros ((cst.nn_hh))
    def run (self):
        r =y # remained vector
        tt =np.array (range (cst.nn_h)) # list of column indices
        ss =set ([]) # extracted column indices
        while True:
            y -pp @ g_hat
            s =np.argmin (abs (pp [:, tt].conj().T @ r))
            ss.add(s)
            pp_ss =pp [:, list(ss)]
            pp_ss_ps_inv =np.linalg.inv (pp_ss.conj().T @ pp_ss) @ pp_ss.conj().T
            r =y -pp_ss.conj().T @ pp_ss_ps_inv @ y
            if np.linalg.norm(r, 2) <self.eG:
                break

class Ll_ss:
    def __init__(self, pp, y):
        self.pp = pp
        self.y = y
        self.g_hat = np.zeros ((cst.nn_hh))
    def run (self):
        g_hat =numpy.linalg.inv (self.pp @ self.pp.conj().T).inv @ pp.conj().T @ self.y

