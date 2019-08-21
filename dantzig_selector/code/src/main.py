#! /usr/bin/env python3

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time

import constants as cst
import classes as cls
import functions as fct

m_g = int ((cst.num_try_gamma - 1) / 2) # hold
arr_gamma = (np.sqrt (np.log2 (cst.nn_h))
    * np.array (
        [2 ** x for x in reversed (range (m_g, -m_g-1, -1))]))

m_s = int ((cst.num_try_sigma - 1) / 2) # hold
arr_sigma = (np.array (
        [np.sqrt(2) ** x for x in reversed (range (m_s, -m_s-1, -1))]))

lst_legend = []
lst_legend.append ("DS, theoretical")
for gamma in arr_gamma:
    lst_legend.append ("DS, γ = " + '%.2f' % gamma)

ddss_bound = (np.log2 (
    4 * np.sqrt (2 * cst.ll * cst.nn_hh * np.log2 (cst.nn_hh))))
kk = fct.kk()

count_prog = 0
time_start = time.time()
lst_lst_err_abs = [] # each sigma, each method
lst_lst_err_rel = [] # each sigma, each method
for i_sigma in range (cst.num_try_sigma):
    sigma = arr_sigma [i_sigma]
    print ("σ = ", "{0:.2f}".format (sigma), " :", sep = '')
    lst_err_abs = [0] * cst.num_method
    norm_hh = 0
    for _ in range (cst.num_repeat):
        count_prog += 1
        ff_bb =fct.ff_bb()
        ff_rr =fct.ff_rr()
        ww_bb =fct.ww_bb()
        ww_rr =fct.ww_rr()
        hh =fct.hh()
        zz =fct.zz()
        gg = kk.conj().T @ hh @ kk
        norm_hh += np.log2 (np.linalg.norm (hh, ord=2))

        pp = np.kron (
            ff_bb.T @ ff_rr.T @ kk.conj(),
            ww_bb @ ww_rr @ kk)
        g = fct.vectorize (gg)
        z = fct.vectorize (zz)
        y = pp @ g + sigma * z
        i_method = 0

        # Dantzig Selector error bound
        lst_err_abs [i_method] += ddss_bound
        i_method += 1

        # Dantzig Selector
        for gamma in arr_gamma:
            dd_ss = cls.Dd_ss (pp, y, gamma);
            dd_ss.run()
            hh_hat = (kk
                @ fct.inv_vectorize (dd_ss.g_hat, cst.nn_hh, cst.nn_hh)
                @ kk.conj().T)
            lst_err_abs [i_method] += np.log2 (np.linalg.norm (hh_hat -hh, ord=2))
            i_method += 1

        percent_progress = 100 * count_prog / (cst.num_try_sigma * cst.num_repeat)
        print (
            "    experiment ", count_prog, " (",
            "{0:.1f}".format (percent_progress), "%)",
            sep = '', end = '\r')
    lst_err_abs = (np.array (lst_err_abs) / cst.num_repeat).tolist ()
    norm_hh /= cst.num_repeat
    lst_err_rel = (np.array (lst_err_abs) / norm_hh).tolist ()
    lst_lst_err_abs.append (lst_err_abs)
    lst_lst_err_rel.append (lst_err_rel)
    print ("                                ", end = '\r') # clear last line
    print ("    done")
time_stop = time.time ()

print (
    "averaged time elapsed for each experiment: ",
    '%.2f' %
        ((time_stop - time_start) / (60 * cst.num_try_sigma * cst.num_repeat)),
    " (min)")
print (
    "total time elapsed: ",
    '%.2f' % ((time_stop - time_start) / 60),
    " (min)")

arr_x = np.array (np.log2 (arr_sigma))
lst_arr_y_abs = np.array (lst_lst_err_abs).T.tolist () # each method, each sigma
lst_arr_y_abs = [np.array (lst) for lst in lst_arr_y_abs]
lst_arr_y_rel = np.array (lst_lst_err_rel).T.tolist () # each method, each sigma
lst_arr_y_rel = [np.array (lst) for lst in lst_arr_y_rel]

label_x = "Std. of Noise (Log)"
label_y = "Absolute 2-Norm error (Log)"
title = "Channel Estimation Performance, Absolute"
fct.save_table (arr_sigma, lst_arr_y_abs, label_x, label_y, lst_legend, title)
fct.save_plot (arr_sigma, lst_arr_y_abs, label_x, label_y, lst_legend, title)
label_x = "Std. of Noise (Log)"
label_y = "Relative 2-Norm error"
title = "Channel Estimation Performance, Relative"
fct.save_table (arr_sigma, lst_arr_y_rel, label_x, label_y, lst_legend, title)
fct.save_plot (arr_sigma, lst_arr_y_rel, label_x, label_y, lst_legend, title)

