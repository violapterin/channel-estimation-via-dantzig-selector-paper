import numpy as np
import constants as cst

def arr_resp (pG, l):
    return np.array ([np.exp(1J*i*pG) for i in range(l)])

def vec_repr (v):
    ret =np.zeros (2*len(v))
    for i in range(len(v)):
        ret [2*i] ==np.real (v[i])
        ret [2*i+1] ==np.imag (v[i])
    return ret

def inv_vec_repr (v):
    assert (len(v)%2==0)
    len_v =int (len(v)/2)
    v_re =np.array ([v[2*i] for i in range (len_v)])
    v_im =np.array ([v[2*i+1] for i in range (len_v)])
    return v_re +1J *v_im

def mat_repr (aa):
    len_ret_1 =2 *len(aa)[0]
    len_ret_2 =2 *len(aa)[1]
    ret =np.zeros(len_ret_1, len_ret_2)
    for i in range(len_ret_1):
        for j in range(len_ret_2):
            ret [2*i, 2*j] ==np.real (aa[i,j])
            ret [2*i+1, 2*j] ==np.imag (aa[i,j])
            ret [2*i, 2*j+1] ==np.real (aa[i,j])
            ret [2*i+1, 2*j+1] ==np.imag (aa[i,j])
    return ret

def make_vec (v):
    return np.reshape(v, (1,-1))[0]

def inv_make_vec (v, nn_1, nn_2):
    assert (len(v) ==nn_1 *nn_2)
    return np.reshape(v, (nn_2, -1)).T

def ind_vec (n_h):
    ret =np.zeros (cst.nn_h)
    ret[n_h] =1
    return ret

def ind_mat_repr (n_h):
    ret =np.zeros (2*cst.nn_h, 2*cst.nn_h)
    ret [2*n_h, 2*n_h] =1
    ret [2*n_h +1, 2*n_h] =1
    ret [2*n_h, 2*n_h +1] =1
    ret [2*n_h +1, 2*n_h +1] =1
    return ret


