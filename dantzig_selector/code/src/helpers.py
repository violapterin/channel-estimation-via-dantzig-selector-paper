import numpy as np
import constants as cst

def arr_response (pG):
    return np.array ([np.exp (1J *i *pG) for i in range(cst.nn_hh)])

def kk ():
    ret =np.zeros(nn_hh, nn_hh)
    for i in range(nn_h):
        for j in range(nn_h):
            ret[i] [j] =(1 /np.sqrt(nn_h)) *np.exp(2 *np.pi *1J *i *j /nn_hh)
    return ret

def find_repr_vec (v):
    ret =np.zeros (2*len(v))
    for i in range(len(v)):
        ret [2*i] ==np.real (v[i])
        ret [2*i+1] ==np.imag (v[i])
    return ret

def inv_find_repr_vec (v):
    assert (len(v)%2==0)
    len_v =int (len(v)/2)
    v_re =np.array ([v[2*i] for i in range (len_v)])
    v_im =np.array ([v[2*i+1] for i in range (len_v)])
    return v_re +1J *v_im

def find_repr_mat (aa):
    ret =np.zeros (2 *len(aa)[0], 2 *len(aa)[1])
    for i in range (len(aa)[0]):
        for j in range (len(aa)[1]):
            ret [2*i] [2*j] ==np.real (aa[i,j])
            ret [2*i+1] [2*j] ==np.imag (aa[i,j])
            ret [2*i] [2*j+1] ==np.real (aa[i,j])
            ret [2*i+1] [2*j+1] ==np.imag (aa[i,j])
    return ret

def vectorize (v):
    return np.reshape(v, (1,-1))[0]

def inv_vectorize (v, nn_1, nn_2):
    assert (len(v) ==nn_1 *nn_2)
    return np.reshape(v, (nn_2, -1)).T

def indic_vec (n_h):
    ret =np.zeros (cst.nn_h)
    ret[n_h] =1
    return ret

def indic_repr_mat (n_h):
    ret =np.zeros (2*cst.nn_h, 2*cst.nn_h)
    ret [2*n_h] [2*n_h] =1
    ret [2*n_h +1] [2*n_h +1] =1
    return ret


