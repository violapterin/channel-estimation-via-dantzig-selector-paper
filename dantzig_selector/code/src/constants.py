import numpy as np


lG =3
ll =5
d_arr =5
nn_yy =3
nn_rr =5
nn_hh =7
nn_y =nn_yy *nn_yy
nn_h =nn_hh *nn_hh
gG =1

kk =np.zeros(nn_hh, nn_hh)
for i in range(nn_h):
    for j in range(nn_h):
        kk[i,j] =() *np.exp(2 *np.pi *1J *i *j /nn_hh)

