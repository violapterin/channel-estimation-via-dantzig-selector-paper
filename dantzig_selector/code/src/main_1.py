#! /usr/bin/env python3

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time

import constants as cst
import classes as cls
import functions as fct


# # # # # # # # # # # # # # # # 

pow_low_gG = -int (cst.num_try_gG * 2 / 3)
pow_high_gG = cst.num_try_gG + pow_low_gG
list_gG = [
    np.sqrt (np.log (cst.nn_h)) * 2 ** x
    for x in reversed (range (int(pow_high_gG), int(pow_low_gG), -1))]
pow_low_noise = -int (cst.num_scale_noise * 2 / 3)
pow_high_noise = cst.num_scale_noise + pow_low_noise
list_noise = [
    np.sqrt(2) ** x
    for x in reversed (range (int(pow_high_noise), int(pow_low_noise), -1))]

kk = fct.kk()

count = 0
time_start = time.time()
list_arr_e = []
for gG in list_gG:
    print ("gamma = ", "{0:.2f}".format (gG), " :", sep ='')
    arr_e = np.array ([])
    for noise in list_noise:
        print ("    Noise std. = ", "{0:.2f}".format (noise), " :", end='', sep ='')
        count += 1
        e = 0
        for _ in range (cst.num_repeat):
            ff_bb =fct.ff_bb()
            ff_rr =fct.ff_rr()
            ww_bb =fct.ww_bb()
            ww_rr =fct.ww_rr()
            hh =fct.hh()
            zz =fct.zz()
            gg = kk.conj().T @ hh @ kk

            pp = np.kron (
                ff_bb.T @ ff_rr.T @ kk.conj(),
                ww_bb @ ww_rr @ kk)
            g = fct.vectorize (gg)
            z = fct.vectorize (zz)
            y = pp @ g + noise * z

            dd_ss = cls.Dd_ss (pp, y, gG);
            dd_ss.run()
            hh_hat = (kk
                @ fct.inv_vectorize (dd_ss.g_hat, cst.nn_hh, cst.nn_hh)
                @ kk.conj().T)
            e += np.log ((np.linalg.norm (hh_hat -hh, ord='fro') ** 2)
                / (np.linalg.norm (hh, ord='fro') ** 2))
        e /= cst.num_repeat
        arr_e =np.append (arr_e, [e])
        count_progress = 100 * count / (cst.num_scale_noise * cst.num_try_gG)
        print (" done (", "{0:.1f}".format (count_progress), "%)", sep = '')
    list_arr_e.append(arr_e)

time_stop = time.time()

fct.draw (list_noise, list_arr_e, "Std. of Noise", "Log Relative Error", "DS_Performance")

print (
    "averaged time elapsed for each data point: ",
    ((time_stop - time_start)
        / (60 * cst.num_scale_noise * cst.num_repeat * cst.num_try_gG)),
    " (min)")
print (
    "total time elapsed: ",
    (time_stop - time_start) / 60,
    " (min)")
