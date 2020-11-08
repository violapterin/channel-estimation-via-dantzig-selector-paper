import numpy as np
from enum import Enum
import sys

import constants as cst
import functions as fct


class Beamformer:
    def __init__ (s, ver):
        s.ver = ver
        s.pp = np.zeros ((cst.NN_Y (s.ver), cst.NN_H (s.ver)))

    def generate (s):
        kk = fct.get_kk (s.ver)
        ff_bb = fct.pick_mat_bb (s.ver).T
        ff_rr = fct.pick_mat_rr (s.ver).T
        ww_bb = fct.pick_mat_bb (s.ver)
        ww_rr = fct.pick_mat_rr (s.ver)
        pp = np.kron (ff_bb.T @ ff_rr.T @ kk.conj (), ww_bb @ ww_rr @ kk)
        s.pp = pp

class Channel:
    def __init__ (s, pp, s_g, ver):
        s.ver = ver
        s.pp = pp
        s.s_g = s_g
        s.hh = np.zeros ((cst.NN_HH (s.ver), cst.NN_HH (s.ver)), dtype = complex)
        s.y = np.zeros ((cst.NN_Y (s.ver)))

    def generate (s):
        s.hh = fct.pick_hh (s.ver)
        s.norm_hh = np.linalg.norm (s.hh, ord = 'fro')

    def zero (s):
        s.hh = np.zeros ((cst.NN_HH (s.ver), cst.NN_HH (s.ver)), dtype = complex)

    def transmit (s):
        kk = fct.get_kk (s.ver)
        zz = fct.pick_zz (s.ver)
        gg = kk.conj ().T @ s.hh @ kk
        g = fct.vectorize (gg)
        z = fct.vectorize (zz)
        s.y = s.pp @ g + (s.s_g / np.sqrt(2)) * z


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Estimation:
    def __init__ (s, cnt_each, hh, pp, y, s_g, ver):
        s.ver = ver
        s.cnt_each = cnt_each
        s.hh = hh
        s.pp = pp
        s.y = y
        s.s_g = s_g
        s.g_r_h = np.zeros (2 * (cst.NN_H (s.ver)))
        s.hh_h = np.zeros ((cst.NN_HH (s.ver), cst.NN_HH (s.ver)), dtype = complex)
        s.rr = 0

    def zero (s):
        s.g_r_h = np.zeros (2 * (cst.NN_H (s.ver)))
        s.hh_h = np.zeros ((cst.NN_HH (s.ver), cst.NN_HH (s.ver)), dtype = complex)
        s.rr = 0

    def set_g_r_h (s, g_r_h):
        s.g_r_h = g_r_h

    def convert (s):
        g_h = fct.inv_find_rep_vec (s.g_r_h)
        gg_h = fct.inv_vectorize (g_h, cst.NN_HH (s.ver), cst.NN_HH (s.ver))
        s.hh_h = (fct.get_kk (s.ver) @ gg_h @ fct.get_kk (s.ver).conj().T)

        nor_hh = np.linalg.norm (s.hh, ord = 2)
        nor_ee = np.linalg.norm (s.hh - s.hh_h, ord = 2)
        s.rr = nor_ee / nor_hh
        
        if (s.rr > cst.MAX_RR (s.ver)):
            pp_r = find_rep_mat (s.pp)
            y_r = find_rep_vec (s.y)
            pp_r_inv = np.linalg.pinv (pp_r)
            s.g_r_h = pp_r_inv @ y_r
            nor_hh = np.linalg.norm (s.hh, ord = 2)
            nor_ee = np.linalg.norm (s.hh - s.hh_h, ord = 2)
            s.rr = nor_ee / nor_hh

    def find_ddss_theory (s, s_g):
        nor_hh = np.linalg.norm (s.hh, ord = 2)
        nor_ee = (8 * s_g * (cst.LL (s.ver) ** (1/2))
                * (np.log (cst.NN_HH (s.ver)) ** (3/2)))
        s.rr = nor_ee / nor_hh

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Version:
    def __init__ (s, focus, size, ratio, iden):
        s.focus = focus
        s.size = size
        s.ratio = ratio
        s.iden = iden

    def get_num (s):
        focus_dict = {
            Focus.ASSORTED: 0,
            Focus.OOMMPP: 1,
            Focus.DDSS: 2
        }
        size_dict = {
            Size.VERY_SMALL: 0,
            Size.SMALL: 1,
            Size.MEDIUM: 2,
            Size.BIG: 3,
            Size.VERY_BIG: 4
        }
        ratio_dict = {
            Ratio.NARROW: 0,
            Ratio.WIDE: 1
        }

        ret = 0
        ret += (10 * focus_dict [s.focus]
                + 2 * size_dict [s.size] + ratio_dict [s.ratio])
        return str(ret).zfill(2)

class Focus (Enum):
    ASSORTED = 0
    OOMMPP = 1
    DDSS = 2

class Size (Enum):
    VERY_SMALL = 0
    SMALL = 1
    MEDIUM = 2
    BIG = 3
    VERY_BIG = 4

class Ratio (Enum):
    NARROW = 0
    WIDE = 1
