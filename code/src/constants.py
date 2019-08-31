import numpy as np
import classes as cls

def NUM_GRID_PHASE (ver):
    return 16

def LAMBDA_ANT (ver):
    return 3

def DIST_ANT ():
    return 2

def NN_YY (ver):
    switcher = {
        cls.Size.TEST : 3,
        cls.Size.SMALL : 6,
        cls.Size.MEDIUM : 9,
        cls.Size.BIG : 12}
    return switcher [ver.size]

def NN_RR (ver):
    switcher = {
        cls.Size.TEST : 6,
        cls.Size.SMALL : 12,
        cls.Size.MEDIUM : 18,
        cls.Size.BIG : 24}
    return switcher [ver.size]

def NN_HH (ver):
    switcher = {
        cls.Size.TEST : 12,
        cls.Size.SMALL : 24,
        cls.Size.MEDIUM : 36,
        cls.Size.BIG : 48}
    return switcher [ver.size]

def LL (ver):
    switcher = {
        cls.Size.TEST : 1,
        cls.Size.SMALL : 2,
        cls.Size.MEDIUM : 3,
        cls.Size.BIG : 4}
    return switcher [ver.size]

def VALUE_SPACING (ver):
    return np.sqrt (2)

def NUM_SIGMA (ver):
    return 7

def NUM_ETA (ver): # OMP
    if (ver.focus == cls.Focus.OOMMPP):
        return 3
    elif (ver.focus == cls.Focus.ASSORTED):
        return 1
    else: # cls.Focus.DDSS
        return 0

def NUM_GAMMA_DS (ver): # DS
    if (ver.focus == cls.Focus.DDSS):
        return 3
    elif (ver.focus == cls.Focus.ASSORTED):
        return 1
    else: # cls.Focus.OOMMPP
        return 0

def GAMMA_LASSO (ver):
    return np.sqrt (2 * np.log2 (cst.NN_H (ver)))

def GAMMA_DDSS (ver):
    return np.sqrt (2 * np.log2 (cst.NN_H (ver)))

def MAX_ITER_OOMMPP (ver):
    return 4 * NN_HH (ver)

def ETA_OOMMPP_2_NORM (ver):
    return 2 * np.sqrt (3 * cst.NN_YY (ver))

def ETA_OOMMPP_INFTY_NORM (ver):
    return 2 * np.sqrt (np.log (cst.NN_H (ver)))
