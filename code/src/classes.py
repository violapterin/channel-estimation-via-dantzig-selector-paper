import numpy as np
import cvxpy as cp
from sklearn.linear_model import Lasso as lasso
import sys
from enum import Enum

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

class Version:
    def __init__ (self, size, focus):
        self.size = size
        self.focus = focus

class Size (Enum):
    TEST = 1
    SMALL = 2
    MEDIUM = 3
    BIG = 4

class Focus (Enum):
    OOMMPP = 1
    DDSS = 2
    ASSORTED = 3
