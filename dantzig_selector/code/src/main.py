#! /usr/bin/env python3

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
import time

import constants as cst
import algorithms as alg
import functions as fct

m_s = int ((cst.num_try_sigma - 1) / 2) # hold
arr_sigma = (np.array (
        [np.sqrt(2) ** x for x in reversed (range (m_s, -m_s-1, -1))]))
m_g = int ((cst.num_try_gamma - 1) / 2) # hold
arr_gamma_ddss = (np.log2 (cst.nn_hh)
    * np.array ([4 ** x for x in reversed (range (m_g, -m_g-1, -1))]))

lst_legend = []
lst_method = [] # arguments: (estimation, sigma)
# DS
# Least Square
lst_legend.append ("LS")
lst_method.append (lambda x, y: llss (x))
# Lasso
lst_legend.append ("Lasso")
lst_method.append (lambda x, y: lasso (x, np.log2 (cst.nn_hh)))
# Orthogonal Matching Pursuit: fixed iteration number
lst_legend.append ("OMP, $L$ times")
lst_method.append (lambda x, y: oommpp_fixed_times (x, y * cst.nn_ll))
# Orthogonal Matching Pursuit: limited l-2 norm
lst_legend.append ("OMP, $l_2$-norm")
lst_method.append (lambda x, y: oommpp_2_norm (x, y * cst.nn_yy))
# Orthogonal Matching Pursuit: limited l-infinity norm
lst_legend.append ("OMP, $l_\infty$-norm")
lst_method.append (lambda x, y: oommpp_infty_norm (x, y * np.log (cst.nn_hh)))
# Dantzig Selector error bound
lst_legend.append ("DS, theory")
lst_method.append (lambda x, y: ddss_theory (x))
for gamma in arr_gamma_ddss:
    lst_legend += ["DS, γ = " + '%.2f' % gamma]
    lst_method.append (lambda x, y: ddss_complex (x, gamma))

num_method = len (lst_method)
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
        estimation = cls.Estimation (pp, y, hh)
        i_method = 0

        # Dantzig Selector error bound
        lst_err_abs [0] += ddss_bound
        for i in range (num_method):
            lst_method [i] (estimation, sigma)
            lst_err_abs [i_method] += estimation.d



        estimation.oommpp_limited_2_norm (sigma * np.log (cst.nn_hh))
        lst_err_abs [i_method] += estimation.d
        i_method += 1

        lst_err_abs [i_method] += ddss_bound
        i_method += 1

        # Dantzig Selector
        for gamma in arr_gamma_ddss:
            ddss_complex = cls.Ddss_complex (pp, y, gamma)
            ddss_complex.run()
            hh_hat = (kk
                @ fct.inv_vectorize (ddss_complex.g_hat, cst.nn_hh, cst.nn_hh)
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

