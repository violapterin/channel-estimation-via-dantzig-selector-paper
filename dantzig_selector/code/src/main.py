#! /usr/bin/env python3

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time

import constants as cst
import classes as cls
import functions as fct

low_gamma = -int (cst.num_try_gamma * 2 / 3)
high_gamma = cst.num_try_gamma + low_gamma
data_gamma = (np.sqrt (np.log2 (cst.nn_h))
    * np.array (
        [2 ** x
            for x in reversed (range (int(high_gamma), int(low_gamma), -1))]))

low_sigma = -int (cst.num_scale_sigma * 2 / 3)
high_sigma = cst.num_scale_sigma + low_sigma
data_sigma = ((1/8)
    * np.array (
        [np.sqrt(2) ** x
            for x in reversed (range (int(high_sigma), int(low_sigma), -1))]))
data_x = np.array (np.log2 (data_sigma))

list_legend = ["DS bound"]
for gamma in data_gamma:
    list_legend.append ("γ = " + '%.2f' % gamma)

abs_bound = (np.log2 (4 * np.sqrt (2 * cst.ll * cst.nn_hh * np.log2 (cst.nn_hh)))
    * np.array ([1 for _ in range (len (data_sigma))]))
list_data_y_abs =[abs_bound]
rel_bound = (np.log2 (4 * np.sqrt (2 * cst.ll * cst.nn_hh * np.log2 (cst.nn_hh)))
    / np.log2 (cst.nn_hh)
    * np.array ([1 for _ in range (len (data_sigma))]))
list_data_y_rel =[rel_bound]

kk = fct.kk()

count = 0
time_start = time.time()
for gamma in data_gamma:
    print ("γ = ", "{0:.2f}".format (gamma), " :", sep ='')
    data_abs_error = np.array ([])
    data_rel_error = np.array ([])
    for sigma in data_sigma:
        print (
            "    σ = ", "{0:.2f}".format (sigma), " :",
            end = '', sep = '', flush = True)
        count += 1
        error = 0
        norm = 0
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
            y = pp @ g + sigma * z

            dd_ss = cls.Dd_ss (pp, y, gamma);
            dd_ss.run()
            hh_hat = (kk
                @ fct.inv_vectorize (dd_ss.g_hat, cst.nn_hh, cst.nn_hh)
                @ kk.conj().T)
            error += np.log2 (np.linalg.norm (hh_hat -hh, ord=2))
            norm += np.log2 (np.linalg.norm (hh, ord=2))
        error /= cst.num_repeat
        norm /= cst.num_repeat
        data_abs_error = np.append (data_abs_error, np.array([error]))
        data_rel_error = np.append (data_rel_error, np.array([error / norm]))
        count_progress = 100 * count / (cst.num_scale_sigma * cst.num_try_gamma)
        print (" done (", "{0:.1f}".format (count_progress), "%)", sep = '')
    list_data_y_abs.append (data_abs_error)
    list_data_y_rel.append (data_rel_error)
time_stop = time.time()

label_x = "Std. of Noise (Log)"
label_y = "Absolute 2-Norm Error (Log)"
title = "Channel Estimation Performance, Absolute"
fct.save_table (data_sigma, list_data_y_abs, label_x, label_y, list_legend, title)
fct.save_plot (data_sigma, list_data_y_abs, label_x, label_y, list_legend, title)

label_x = "Std. of Noise (Log)"
label_y = "Relative 2-Norm Error"
title = "Channel Estimation Performance, Relative"
fct.save_table (data_sigma, list_data_y_rel, label_x, label_y, list_legend, title)
fct.save_plot (data_sigma, list_data_y_rel, label_x, label_y, list_legend, title)

print (
    "averaged time elapsed for each data point: ",
    '%.2f' %
        ((time_stop - time_start)
            / (60 * cst.num_scale_sigma * cst.num_repeat * cst.num_try_gamma)),
    " (min)")
print (
    "total time elapsed: ",
    '%.2f' % ((time_stop - time_start) / 60),
    " (min)")
