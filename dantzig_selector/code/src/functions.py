import numpy as np
import matplotlib.pyplot as plt
import os

import constants as cst
import classes as cls

def array_response (pG):
    return np.array ([np.exp (1J *i *pG) for i in range (cst.nn_hh)])

def kk (): # DFT matrix
    ret =np.zeros ((cst.nn_hh, cst.nn_hh))
    for i in range (cst.nn_hh):
        for j in range (cst.nn_hh):
            ret [i] [j] =(1 /np.sqrt (cst.nn_h)) *np.exp (2 *np.pi *1J *i *j /cst.nn_hh)
    return ret

def find_repr_vec (v):
    ret =np.zeros ((2*len (v)))
    for i in range (len (v)):
        ret [2*i] ==np.real (v [i])
        ret [2*i+1] ==np.imag (v [i])
    return ret

def inv_find_repr_vec (v):
    assert (len (v)%2==0)
    len_v =int (len (v)/2)
    v_re =np.array ([v [2*i] for i in range (len_v)])
    v_im =np.array ([v [2*i+1] for i in range (len_v)])
    return v_re +1J *v_im

def find_repr_mat (aa):
    ret =np.zeros ((2 *len (aa)[0], 2 *len (aa)[1]))
    for i in range (len (aa)[0]):
        for j in range (len (aa)[1]):
            ret [2*i] [2*j] ==np.real (aa [i,j])
            ret [2*i+1] [2*j] ==np.imag (aa [i,j])
            ret [2*i] [2*j+1] ==np.real (aa [i,j])
            ret [2*i+1] [2*j+1] ==np.imag (aa [i,j])
    return ret

def vectorize (v):
    return np.reshape (v, (1,-1))[0]

def inv_vectorize (v, nn_1, nn_2):
    assert (len (v) ==nn_1 *nn_2)
    return np.reshape (v, (nn_2, -1)).T

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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def draw (arr_x, list_dat, label_x, label_y, title):
    plt.close ("all")
    fig =plt.figure ()
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
    size_marker = 7 # fixed marker size
    width_line = 3 # fixed line width

    size_line =len (arr_x)
    for i in range (len (list_dat)):
        dat =list_dat[i]
        assert (size_line == len (dat.arr))
        plt.plot (arr_x,
                  dat.arr,
                  markersize =size_marker,
                  linewidth =width_line,
                  linestyle =list_style [int (i % num_style)],
                  color =list_color [int (i % num_color)],
                  marker =list_marker [int (i % num_marker)],
                  )
    # TODO: legend

    path_fig_out =(os.path.abspath (os.path.join (os.getcwd (), os.path.pardir))
                  +"/plt/" +title +".png")
    if os.path.isfile(path_fig_out):
        os.system("rm "+path_fig_out)
    plt.savefig (path_fig_out, bbox_inches ="tight")





