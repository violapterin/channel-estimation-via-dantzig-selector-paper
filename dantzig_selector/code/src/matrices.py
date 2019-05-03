import numpy as np
import constants as cst

class Hh:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_h, cst.nn_h))
    
    def generate (self):
        self.val = np.zeros ((cst.nn_h, cst.nn_h))
        for l in range (cst.ll_path):
            aG_l = np.random.normal (cst.amp_mean_hh, cst.amp_std_hh)
            pG_l = (dist_antenna /cst.lG) * np.sin (np.random.uniform (0, 2 * np.pi))
            tG_l = (dist_antenna /cst.lG) * np.sin (np.random.uniform (0, 2 * np.pi))
            self.val += aG_l * (array_response (pG_l)) @ array_response (pG_l).H

class Zz:
    def __init__ (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_yy))
    
    def generate (self):
        self.val = np.zeros ((cst.nn_yy, cst.nn_yy))
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
        
