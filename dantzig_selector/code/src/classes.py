import numpy as np
import cvxpy as cp

import constants as cst
import functions as fct

class Hh:
    def __init__ (self):
        self.zero ()

    def zero (self):
        self.val = np.zeros ((cst.nn_h, cst.nn_h), dtype=complex)
    
    def generate (self):
        self.val = np.zeros ((cst.nn_h, cst.nn_h), dtype=complex)
        for l in range (cst.ll_path):
            aG_l = np.random.normal (cst.amp_mean_hh, cst.amp_std_hh)
            pG_l = (cst.dist_antenna /cst.lG) * np.sin (np.random.uniform (0, 2 * np.pi))
            tG_l = (cst.dist_antenna /cst.lG) * np.sin (np.random.uniform (0, 2 * np.pi))
            self.val += aG_l * (fct.array_response (pG_l)) @ fct.array_response (pG_l).conj().T

class Zz:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_yy), dtype=complex)
    
    def generate (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_yy), dtype=complex)
        for i in range (cst.nn_yy):
            for j in range (cst.nn_yy):
                self.val [i] [j] = (
                    np.random.normal (0, 1)
                    + 1J * np.random.normal (0, 1))

class Ff_bb:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_rr, cst.nn_yy))

    def generate (self):
        self.val = np.zeros ((cst.nn_rr, cst.nn_yy))
        for i in range (cst.nn_rr):
            for j in range (cst.nn_yy):
                self.val [i] [j] = (
                    np.random.normal (cst.part_mean_bb, cst.part_std_bb)
                    + 1J * np.random.normal (cst.part_mean_bb, cst.part_std_bb))
        self.val /= power_constr

class Ff_rr:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_hh, self.cst.nn_rr))

    def generate (self):
        for i in range (cst.nn_hh):
            for j in range (cst.nn_rr):
                idx_phase = np.random.randint (cst.num_grid_phase)
                self.val[i][j] = np.exp (2 * np.pi * 1J * idx_phase /num_grid_phase)
        
class Ww_bb:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_rr))

    def generate (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_rr))
        for i in range (cst.nn_yy):
            for j in range (cst.nn_rr):
                self.val [i] [j] = (
                        np.random.normal (cst.part_mean_bb, cst.part_std_bb)
                        + 1J * np.random.normal (cst.part_mean_bb, cst.part_std_bb))
        self.val /= power_constr

class Ww_rr:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_rr, self.cst.nn_hh))

    def generate (self):
        for i in range (cst.nn_rr):
            for j in range (cst.nn_hh):
                idx_phase = np.random.randint (cst.num_grid_phase)
                self.val[i][j] = np.exp (2 * np.pi * 1J * idx_phase /num_grid_phase)
        

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

