import numpy as np
import cvxpy as cp
import constants as cst

class Ddss(cst):
    def __init__(self):
        self.pp = np.ones ((cst.nn_y, cst.nn_h))
        self.y = np.zeros ((cst.nn_y))
        self.gG = float('inf')
        self.g_h = np.zeros ((cst.nn_h))
        self.hh_h = np.zeros ((cst.nn_hh, cst.nn_hh))

    def set(self, pp, y, gG):
        self.pp = pp
        self.y = y
        self.gG = gG

    def run(self):

    t =np.concatenate ([np.zeros ((2*cst.nn_h)), np.ones ((cst.nn_h))])

for i in range(cst.nn_h):
    A.append (np.concatenate ([ind_mat_repr (i), +np.zeros ((cst.nn_h))]))
    b.append (np.zeros ((2*cst.nn_h)))
    c.append (np.concatenate ([np.zeros ((2*cst.nn_h)), ind_vec (i)]))
    d.append (0)

for i in range(cst.nn_h, 2*cst.nn_h):
    A.append (np.concatenate ([ind_mat_repr (i -cst.nn_h)
                                   @ mat_repr (self.pp) .H
                                   @ mat_repr (self.pp),
                               np.zeros ((cst.nn_h))]), axis=1)
    A.append (ind_mat_repr (i -cst.nn_h)
                  @ mat_repr (self.pp) .H
                  @ mat_repr (self.y))
    c.append (np.concatenate (np.zeros ((3*cst.nn_h))))
    d.append (cst.gG)

x = cp.Variable (3*cst.nn_h)
constr = [cp.SOC (c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range (2*cst.nn_h)]
prob = cp.Problem (cp.Minimize(f.T@x), constr)
prob.solve()

x_hat =x.value
g_repr_hat =x_hat [0:2*cst.nn_h -1]
g_hat =inv_vec_repr (g_repr_hat)
gg_hat =inv_vec (g_hat)
hh_hat =cst.kk @ gg_hat @ cst.kk.H


# Print result.
print("The optimal value is", prob.value)
print("A solution x is")
print(x.value)
for i in range(m):
    print("SOC constraint %i dual variable solution" % i)
    print(soc_constraints[i].dual_value)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#class Oommpp:
#    def __init__(self):
#        pass

