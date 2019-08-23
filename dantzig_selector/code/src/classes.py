import numpy as np
import cvxpy as cp

import constants as cst
import functions as fct
import sys


class Ddss:
    def __init__(self, pp, y, gamma):
        self.pp = pp
        self.y = y
        self.gamma = gamma
        self.g_hat = np.zeros ((cst.nn_h))

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
            d.append (self.gamma)
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

class Oommpp:
    def __init__(self, pp, y, epsilon):
        self.pp = pp
        self.y = y
        self.epsilon = epsilon
        self.g_hat = np.zeros ((cst.nn_h), dtype = complex)
    def run (self):
        r = self.y # remained vector
        tt = range (cst.nn_h) # list of column indices
        ss = [] # extracted column indices
        count_iter = 0
        pp_ss_inv = np.zeros ((cst.nn_h, cst.nn_y))
        while True:
            lst_match = [abs (self.pp [:, i].conj().T @ r) for i in tt]
            s = np.argmax (lst_match)
            ss.append (s)
            ss = sorted (ss)
            pp_ss = self.pp [:, ss]
            if (abs (np.linalg.det (pp_ss.conj().T @ pp_ss)) <= 1.0e-08):
                sys.exit("Error: Singular matrix!")
            pp_ss_inv = np.linalg.inv (pp_ss.conj().T @ pp_ss) @ pp_ss.conj().T
            r = self.y - pp_ss @ pp_ss_inv @ self.y
            count_iter += 1
            if ((np.linalg.norm (r, 2) <self.epsilon)
                    or (count_iter >= cst.max_iter_oommpp)):
                break
        g_hat_ss = pp_ss_inv @ self.y
        for i in range (len(ss)):
            self.g_hat [ss[i]] = g_hat_ss[i] 

class Llss:
    def __init__(self, pp, y):
        self.pp = pp
        self.y = y
        self.g_hat = np.zeros ((cst.nn_h), dtype = complex)
    def run (self):
        self.g_hat = (self.pp.conj().T @
            np.linalg.inv (self.pp @ self.pp.conj().T) @ self.y)

