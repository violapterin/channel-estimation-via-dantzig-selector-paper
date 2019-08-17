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

def save_plot (data_x, list_data_y, label_x, label_y, list_legend, title):
    plt.close ("all")
    plt.title (title, fontsize = 15)
    plt.xlabel (label_x, fontsize = 12)
    plt.ylabel (label_y, fontsize = 12)

    num_style = 4
    list_style = ['-', '--', '-.', ':']
        # from the thinner to thicker
    num_color = 5
    list_color = ['r', 'g', 'c', 'b', 'k']
        # red, green, cyan, blue, black
    num_marker = 6
    list_marker = ['v', '^', 'o', 's', '*', 'D']
        # triangle down, triangle up, circle, square, star, diamond
    size_marker = 7
    width_line = 3

    assert (len (list_data_y) == len (list_legend))
    size_line =len (data_x)
    i = 0
    for data_y in list_data_y:
        assert (size_line == len (data_y))
        plt.plot (
            data_x,
            data_y,
            markersize = size_marker,
            linewidth = width_line,
            linestyle = list_style [int (i % num_style)],
            color = list_color [int (i % num_color)],
            marker = list_marker [int (i % num_marker)],
            label = list_legend [i])
        i += 1
    plt.legend (
        bbox_to_anchor = (1.05, 1),
        loc = 'upper left',
        borderaxespad = 0.)

    path_plot_out =(
        os.path.abspath (os.path.join (os.getcwd (), os.path.pardir))
        + "/plt/" + title.replace (" ", "_") +".png")
    if os.path.isfile(path_plot_out):
        os.system("rm "+path_plot_out)
    plt.savefig (path_plot_out, bbox_inches = "tight")

def save_table (data_x, list_data_y, label_x, label_y, list_legend, title):

    path_table_out =(
        os.path.abspath (os.path.join (os.getcwd (), os.path.pardir))
        + "/dat/" + title.replace (" ", "_") + ".txt")

    with open(path_table_out, 'w') as the_file:
        the_file.write (label_x + '\t')
        the_file.write ('\t'.join (map (str, data_x)) + '\n')
        i = 0
        for data_y in list_data_y:
            the_file.write (
                list_legend [i] + '\t'
                    + '\t'.join (map (str, list_data_y[i])) + '\n')
            i += 1




