#! /usr/bin/env python3

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time

import constants as cst
import functions as fct
import algorithms as alg

arr_sigma = np.array (
    [2 ** x for x in list (
        (np.array(range (cst.num_try_sigma))
            - (cst.num_try_sigma - 1) / 2) / 2)])
arr_gamma_0_ddss = (np.sqrt (2 * np.log2 (cst.nn_h))
    * np.array (
        [2 ** x for x in list (
            (np.array(range (cst.num_try_gamma))
                - (cst.num_try_gamma - 1) / 2) / 2)]))
lst_legend = []
lst_method = [] # arguments: (estimation, sigma)
# Least Square
lst_legend.append ("LS")
lst_method.append (
    lambda x, y: alg.llss (x))
# Lasso
lst_legend.append ("Lasso")
lst_method.append (
    lambda x, y: alg.lasso (x, np.sqrt (2 * np.log2 (cst.nn_h))))
# Orthogonal Matching Pursuit: fixed iteration number
lst_legend.append ("OMP, $L$ times")
lst_method.append (
    lambda x, y: alg.oommpp_fixed_times (x, y * cst.ll))
# Orthogonal Matching Pursuit: limited l-2 norm
lst_legend.append ("OMP, $l_2$-norm")
lst_method.append (
    lambda x, y: alg.oommpp_2_norm (x, 2 * y * np.sqrt (3 * cst.nn_yy)))
# Orthogonal Matching Pursuit: limited l-infinity norm
lst_legend.append ("OMP, $l_\infty$-norm")
lst_method.append (
    lambda x, y: alg.oommpp_infty_norm (x, y * 2 * np.sqrt (2 * np.log (cst.nn_hh))))
# Dantzig Selector error bound
lst_legend.append ("DS, theory")
lst_method.append (
    lambda x, y: alg.ddss_theory (x))
for gamma_0 in arr_gamma_0_ddss:
    lst_legend.append ("DS, $\gamma$ = " + '%.2f' % gamma_0 + "$\sigma$")
    lst_method.append (
        lambda x, y, g = gamma_0: alg.ddss_complex (x, np.sqrt (2) * y * g))
assert (len (lst_method) == len (lst_legend))
num_method = len (lst_method)

count_prog = 0
time_start = time.time()
lst_lst_err_abs = [] # each sigma, each method
lst_lst_err_rel = [] # each sigma, each method
for i_sigma in range (cst.num_try_sigma):
    sigma = arr_sigma [i_sigma]
    print ("σ = " + '%.2f' % sigma + " :", sep = '')
    lst_err_abs = [0] * num_method
    norm_hh = 0
    for _ in range (cst.num_repeat):
        count_prog += 1
        ff_bb = fct.ff_bb()
        ff_rr = fct.ff_rr()
        ww_bb = fct.ww_bb()
        ww_rr = fct.ww_rr()
        hh = fct.hh()
        zz = fct.zz()
        gg = fct.kk().conj().T @ hh @ fct.kk()
        norm_hh += np.log2 (np.linalg.norm (hh, ord=2))

        pp = np.kron (
            ff_bb.T @ ff_rr.T @ fct.kk().conj(),
            ww_bb @ ww_rr @ fct.kk())
        g = fct.vectorize (gg)
        z = fct.vectorize (zz)
        y = pp @ g + sigma * z
        est = alg.Estimation (pp, y, hh)

        for i in range (num_method):
            lst_method [i] (est, sigma)
            lst_err_abs [i] += np.log2 (est.d)

        percent_progress = 100 * count_prog / (cst.num_try_sigma * cst.num_repeat)
        print (
            "    experiment " + count_prog + " ("
            + '%.1f' % percent_progress + "%)",
            sep = '', end = '\r')
    lst_err_abs = list ((np.array (lst_err_abs) / cst.num_repeat))
    norm_hh /= cst.num_repeat
    lst_err_rel = list ((np.array (lst_err_abs) / norm_hh))
    lst_lst_err_abs.append (lst_err_abs)
    lst_lst_err_rel.append (lst_err_rel)
    print ("                                ", end = '\r') # clear last line
    print ("    done")
time_stop = time.time ()

print (
    "averaged time elapsed for each experiment: "
    + '%.2f' %
        ((time_stop - time_start) / (60 * cst.num_try_sigma * cst.num_repeat))
    + " (min)")
print (
    "total time elapsed: ",
    + '%.2f' % ((time_stop - time_start) / 60)
    + " (min)")

arr_x = np.array (np.log2 (arr_sigma))
lst_arr_y_abs = list (np.array (lst_lst_err_abs).T) # each method, each sigma
lst_arr_y_abs = [np.array (lst) for lst in lst_arr_y_abs]
lst_arr_y_rel = list (np.array (lst_lst_err_rel).T) # each method, each sigma
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

