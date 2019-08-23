#! /usr/bin/env python3

import matplotlib.pyplot as plt # plotting functions
import numpy as np
import constants as cst
import classes as cls
import os # getcwd
#import cvxpy as cp

import functions as fct

a = [3,5,0]
a.append(1)
sorted (a)
print (a)
a = sorted (a)
print (a)
quit ()

pp = np.array ([[7,8,0], [4,5,2*1J]])
g = np.array ([5,9,-3])
z = np.array ([0.1,-0.2])
y = pp @ g + z
ss = set ([0,2])
#lst_ip = [pp [:, i].conj().T @ z for i in list(ss)]
#lst_ip = [np.array (pp [:, i]).conj().T @ z for i in list(ss)]
print (pp [:, 2].conj())
quit()

pp_ss = pp [:, list(ss)]
pp_ss_inv = pp_ss.conj().T @ np.linalg.inv (pp_ss @ pp_ss.conj().T)
g_ss = pp_ss_inv @ y
g = np.zeros (3)
for i in range (len(ss)):
    g [sorted(ss)[i]] = g_ss[i] 

print (g_ss)
print (g)

#ll_ss = Ll_ss (pp, y)
#ll_ss. run

# dat_x = np.array ([2,4,6,8])
# dat_y1 =np.zeros (4)
# dat_y2 =np.zeros (4)
# for i in range(4):
#     dat_y1 [i] = 3 * dat_x [i] + 0.1 * np.random.normal (0, 1)
#     dat_y2 [i] = 5 * dat_x [i] + 0.1 * np.random.normal (0, 1)
# list_dat_y =[dat_y1, dat_y2]
# list_legend =["a", "b"]
# fct.save_plot (dat_x, list_dat_y, "foo", "bar", list_legend, "Quux")
# fct.save_table (dat_x, list_dat_y, "foo", "bar", list_legend, "Quux")

#tt =np.array (range(3))
#pp =np. array ([[4,5,6,3], [0,-2,7,4]])
#g =np. array ([4,4,5,5])
#y =pp @ g
#ss =set ([0,3])
#r =y
#t =np.argmin (abs (pp [:, tt].conj().T @ r))



