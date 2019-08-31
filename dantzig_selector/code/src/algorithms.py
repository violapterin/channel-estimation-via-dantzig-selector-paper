import numpy as np
import cvxpy as cp
from sklearn.linear_model import Lasso as lasso
import sys

import constants as cst
import functions as fct

class Estimation:
    def __init__ (self, pp, y, hh):
        self.pp = pp
        self.y = y
        self.hh = hh
        self.g_hat = np.zeros ((cst.nn_h), dtype = complex)
        self.hh_hat = np.zeros ((cst.nn_hh, cst.nn_hh), dtype = complex)
        self.d = 0

    def refresh (self):
        self.g_hat = np.zeros ((cst.nn_h), dtype = complex)
        self.hh_hat = np.zeros ((cst.nn_hh, cst.nn_hh), dtype = complex)
        self.d = 0

    def convert (self):
        self.hh_hat = (fct.kk()
            @ fct.inv_vectorize (self.g_hat, cst.nn_hh, cst.nn_hh)
            @ fct.kk().conj().T)
        self.d = np.linalg.norm (self.hh_hat - self.hh, ord=2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def llss (est):
    est.refresh ()
    try:
        pp_inv = np.linalg.pinv (est.pp)
    except np.linalg.LinAlgError:
        print ("Least square fails due to singularity!", flush = True)
        return
    est.g_hat = pp_inv @ est.y
    est.convert ()

def lasso (est, gamma):
    est.refresh ()
    g = cp.Variable (cst.nn_h, complex = True)
    prob = cp.Problem (
        cp.Minimize (cp.norm (est.pp @ g - est.y, 2)),
        [cp.norm (g, 1) <= gamma])
    try:
        prob.solve ()
    except cp.error.SolverError:
        print ("Lasso fails to solve the program!", flush = True)
        est.g_hat = np.linalg.pinv (est.pp) @ est.y
        return
    est.g_hat = g.value
    est.convert ()

def oommpp_fixed_times (est, times):
    est.refresh ()
    r = est.y # remainder
    tt = range (cst.nn_h) # list of column indices
    ss = [] # extracted column indices
    count_iter = 0
    pp_ss_inv = np.zeros ((cst.nn_h, cst.nn_y))
    while True:
        count_iter += 1
        lst_match = [abs (est.pp [:, i].conj().T @ r) for i in tt]
        s = np.argmax (lst_match)
        ss.append (s)
        ss = list (sorted (set (ss)))
        pp_ss = est.pp [:, ss]
        try:
            pp_ss_inv = np.linalg.pinv (pp_ss)
        except np.linalg.LinAlgError:
            print ("Orthogonal mathcing pursuit fails due to singularity!", flush = True)
            return
        r = est.y - pp_ss @ pp_ss_inv @ est.y
        if (count_iter >= times):
            break
    g_hat_ss = pp_ss_inv @ est.y
    for i in range (len(ss)):
        est.g_hat [ss[i]] = g_hat_ss[i] 
    est.convert ()

def oommpp_2_norm (est, alpha):
    est.refresh ()
    r = est.y # remained vector
    tt = range (cst.nn_h) # list of column indices
    ss = [] # extracted column indices
    count_iter = 0
    pp_ss_inv = np.zeros ((cst.nn_h, cst.nn_y))
    while True:
        count_iter += 1
        lst_match = [abs (est.pp [:, i].conj().T @ r) for i in tt]
        s = np.argmax (lst_match)
        ss.append (s)
        ss = list (sorted (set (ss)))
        pp_ss = est.pp [:, ss]
        try:
            pp_ss_inv = np.linalg.pinv (pp_ss)
        except np.linalg.LinAlgError:
            print ("Orthogonal mathcing pursuit fails due to singularity!", flush = True)
            return
        r = est.y - pp_ss @ pp_ss_inv @ est.y
        if (np.linalg.norm (r, 2) <= alpha):
            break
        if (count_iter >= cst.max_iter_oommpp):
            break
    g_hat_ss = pp_ss_inv @ est.y
    for i in range (len(ss)):
        est.g_hat [ss[i]] = g_hat_ss[i] 
    est.convert ()

def oommpp_infty_norm (est, alpha):
    est.refresh ()
    r = est.y # remained vector
    tt = range (cst.nn_h) # list of column indices
    ss = [] # extracted column indices
    count_iter = 0
    pp_ss_inv = np.zeros ((cst.nn_h, cst.nn_y))
    while True:
        count_iter += 1
        lst_match = [abs (est.pp [:, i].conj().T @ r) for i in tt]
        s = np.argmax (lst_match)
        ss.append (s)
        ss = list (sorted (set (ss)))
        pp_ss = est.pp [:, ss]
        try:
            pp_ss_inv = np.linalg.pinv (pp_ss)
        except np.linalg.LinAlgError:
            print ("Orthogonal mathcing pursuit fails due to singularity!", flush = True)
            return
        r = est.y - pp_ss @ pp_ss_inv @ est.y
        if (np.linalg.norm (pp_ss.conj().T @ r, 2) <= alpha):
            break
        if (count_iter >= cst.max_iter_oommpp):
            break
    g_hat_ss = pp_ss_inv @ est.y
    for i in range (len(ss)):
        est.g_hat [ss[i]] = g_hat_ss[i] 
    est.convert ()

def ddss_theory (est):
    est.refresh ()
    est.d = (np.log2 (4 * np.sqrt (2 * cst.ll * np.log2 (cst.nn_hh))))

def ddss_complex (est, gamma):
    est.refresh ()
    g = cp.Variable (cst.nn_h, complex = True)
    prob = cp.Problem (
        cp.Minimize (cp.norm (g, 1)),
        [cp.norm (est.pp.conj().T @ (est.pp @ g - est.y), "inf")
            <= gamma])
    try:
        prob.solve ()
    except cp.error.SolverError:
        print ("Complex DS fails to solve the program!", flush = True)
        est.g_hat = (np.linalg.pinv (est.pp) @ est.y)
        return
    est.g_hat = g.value
    est.convert ()

# Deprecated!
def ddss_real (est, gamma):
    est.refresh ()
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
                        @ fct.find_repr_mat (est.pp.conj().T)
                        @ fct.find_repr_mat (est.pp),
                    np.zeros ((2 * cst.nn_h, cst.nn_h))],
                axis= 1))
        b.append (fct.indication_repr_mat (i)
            @ fct.find_repr_mat (est.pp.conj().T)
            @ fct.find_repr_vec (est.y))
        c.append (np.zeros ((3 * cst.nn_h)))
        d.append (gamma)
    t = np.concatenate (
        [np.zeros ((2 * cst.nn_h)),
            np.ones ((cst.nn_h))])
    x = cp.Variable (3 * cst.nn_h) # real

    constr = [cp.SOC (c[i].T @ x + d[i], aa[i] @ x + b[i]) for i in range (2 * cst.nn_h)]
    prob = cp.Problem (cp.Minimize (t.T@x), constr)
    try:
        prob.solve ()
    except cp.error.SolverError:
        print ("Real DS fails to solve the program!", flush = True)
        est.g_hat = (np.linalg.pinv (est.pp) @ est.y)
        return
    est.g_hat = g.value

    x_hat = x.value
    g_repr_hat = x_hat [0 : 2 * cst.nn_h]
    est.g_hat = fct.inv_find_repr_vec (g_repr_hat)

    est.convert ()

