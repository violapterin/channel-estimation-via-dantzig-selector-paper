import numpy as np
import constants as cst

def arr_resp(pG, l)
    return np.array ([np.exp(1J*i*pG) for i in range(l)])

def repr (v)
    ret =np.zeros (2*len(v))
    for i in range(len(v))
        ret [2*i] ==np.real (v[i])
        ret [2*i+1] ==np.imag (v[i])
    return ret

def inv_repr (v)
    assert (len(v)%2==0)
    v_re =np.array ([v[2*i] for i in range (int (len(v)/2))])
    v_im =np.array ([v[2*i+1] for i in range (int (len(v)/2))])
    return v_re +1J *v_im

def vec (v)
    return np.reshape(v, (1,-1))[0]

def inv_vec (v, nn_1, nn_2)
    assert (len(v) ==nn_1 *nn_2)
    return np.reshape(v, (nn_2, -1)).T

def ind_vec (n_h)
    ret =np.zeros (cst.nn_h)
    ret[n_h] =1
    return ret

def ind_mat_repr (n_h_1, n_h_2)
    ret =np.zeros (2*cst.nn_h, 2*cst.nn_h)
    ret [2*n_h_1, 2*n_h_2] =1
    ret [2*n_h_1 +1, 2*n_h_2] =1
    ret [2*n_h_1, 2*n_h_2 +1] =1
    ret [2*n_h_1 +1, 2*n_h_2 +1] =1
    return ret


