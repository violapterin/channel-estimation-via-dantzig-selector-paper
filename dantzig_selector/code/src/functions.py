import numpy as np
import matplotlib.pyplot as plt
import os

import constants as cst
import classes as cls

def mat_complex_normal (nn_1, nn_2):
    return ((np.random.normal (0, 1, (nn_1, nn_2))
        +1J * np.random.normal (0, 1, (nn_1, nn_2))))

def mat_uniform_phase (nn_1, nn_2):
    return (
        np.exp (
            2 * np.pi * 1J * np.random.randint (
                cst.num_grid_phase,
                size=(nn_1, nn_2))
            / cst.num_grid_phase))

def zz ():
    return mat_complex_normal (cst.nn_yy, cst.nn_yy)

def ff_bb ():
    return (mat_complex_normal (cst.nn_rr, cst.nn_yy)
    / np.sqrt (cst.nn_yy))

def ww_bb ():
    return (mat_complex_normal (cst.nn_yy, cst.nn_rr)
    / np.sqrt (cst.nn_yy))

def ff_rr ():
    return (mat_uniform_phase (cst.nn_hh, cst.nn_rr)
    / np.sqrt (cst.nn_rr))

def ww_rr ():
    return (mat_uniform_phase (cst.nn_rr, cst.nn_hh)
    / np.sqrt (cst.nn_rr))

def hh():
    ret =np.zeros ((cst.nn_hh, cst.nn_hh), dtype=complex)
    for _ in range (cst.ll):
        aG = (np.random.normal (0, cst.nn_hh / cst.ll)
            + 1J * np.random.normal (0, cst.nn_hh / cst.ll))
        pG = (2 * np.pi * (cst.d_ant /cst.lG_ant)
            * np.sin (np.random.uniform (0, 2 * np.pi)))
        tG = (2 * np.pi * (cst.d_ant /cst.lG_ant)
            * np.sin (np.random.uniform (0, 2 * np.pi)))
        ret += aG * arr_resp (pG) @ arr_resp (tG).conj().T
    return ret


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def kk (): # DFT matrix
    ret =np.zeros ((cst.nn_hh, cst.nn_hh), dtype=complex)
    for i in range (cst.nn_hh):
        for j in range (cst.nn_hh):
            ret [i] [j] =(1 /np.sqrt (cst.nn_hh)) *np.exp (2 *np.pi *1J *i *j /cst.nn_hh)
    return ret

def arr_resp (t):
    return ((1 / np.sqrt (cst.nn_hh))
        * np.array ([np.exp (1J * i * t) for i in range (cst.nn_hh)]))

def find_repr_vec (v):
    ret =np.zeros ((2*len (v)))
    for i in range (len (v)):
        ret [2*i] =np.real (v [i])
        ret [2*i+1] =np.imag (v [i])
    return ret

def inv_find_repr_vec (v):
    assert (len (v)%2 == 0)
    len_v =int (len (v)/2)
    v_re =np.array ([v [2*i] for i in range (len_v)])
    v_im =np.array ([v [2*i+1] for i in range (len_v)])
    return v_re +1J *v_im

def find_repr_mat (aa):
    ret =np.zeros ((2 *(aa.shape[0]), 2 *(aa.shape[1])))
    for i in range (aa.shape[0]):
        for j in range (aa.shape[1]):
            ret [2*i] [2*j] =np.real (aa [i,j])
            ret [2*i+1] [2*j] =np.imag (aa [i,j])
            ret [2*i] [2*j+1] =-np.imag (aa [i,j])
            ret [2*i+1] [2*j+1] =np.real (aa [i,j])
    return ret

def vectorize (v):
    return np.reshape (v, (1,-1))[0]

def inv_vectorize (v, nn_1, nn_2):
    assert (len (v) == nn_1 * nn_2)
    return np.reshape (v, (nn_1, -1))

def indication_vec (n_h):
    ret =np.zeros ((cst.nn_h))
    ret [n_h] =1
    return ret

def indication_repr_mat (n_h):
    ret =np.zeros ((2*cst.nn_h, 2*cst.nn_h))
    ret [2*n_h] [2*n_h] =1
    ret [2*n_h +1] [2*n_h +1] =1
    return ret

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def save_plot (arr_x, lst_arr_y, label_x, label_y, lst_legend, title):
    plt.close ("all")
    plt.title (title, fontsize = 15)
    plt.xlabel (label_x, fontsize = 12)
    plt.ylabel (label_y, fontsize = 12)

    num_style = 4
    lst_style = ['-', '--', '-.', ':']
        # from the thinner to thicker
    num_color = 5
    lst_color = ['r', 'g', 'c', 'b', 'k']
        # red, green, cyan, blue, black
    num_marker = 6
    lst_marker = ['v', '^', 'o', 's', '*', 'D']
        # triangle down, triangle up, circle, square, star, diamond
    size_marker = 7
    width_line = 3

    assert (len (lst_arr_y) == len (lst_legend))
    for i_method in range (len (lst_arr_y)):
        arr_y = lst_arr_y [i_method]
        assert (len (arr_x) == len (arr_y))
        plt.plot (
            arr_x,
            arr_y,
            markersize = size_marker,
            linewidth = width_line,
            linestyle = lst_style [int (i_method % num_style)],
            color = lst_color [int (i_method % num_color)],
            marker = lst_marker [int (i_method % num_marker)],
            label = lst_legend [i_method])
    plt.legend (
        bbox_to_anchor = (1.05, 1),
        loc = 'upper left',
        borderaxespad = 0.)

    path_plot_out =(
        os.path.abspath (os.path.join (os.getcwd (), os.path.pardir))
            + "/plt/" + title.replace (" ", "_") + ".png")
    if os.path.isfile(path_plot_out):
        os.system ("rm " + path_plot_out)
    plt.savefig (path_plot_out, bbox_inches = "tight")

def save_table (arr_x, lst_arr_y, label_x, label_y, lst_legend, title):
    path_table_out =(
        os.path.abspath (os.path.join (os.getcwd (), os.path.pardir))
        + "/dat/" + title.replace (" ", "_") + ".txt")

    with open (path_table_out, 'w') as the_file:
        the_file.write (label_x + '\t')
        the_file.write ('\t'.join (map (str, arr_x)) + '\n')
        assert (len (lst_arr_y) == len (lst_legend))
        for i in range (len (lst_arr_y)):
            the_file.write (
                lst_legend [i] + '\t'
                    + '\t'.join (map (str, lst_arr_y [i])) + '\n')


