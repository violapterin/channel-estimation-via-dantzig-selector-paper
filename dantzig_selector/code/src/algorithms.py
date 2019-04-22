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

t =np.zeros ((2*cst.nn_y)) +np.zeros ((cst.nn_h))

for i in range(nn_h):
    A.append ()
    b.append ()
    c.append ()
    d.append ()

for i in range(nn_h):
    A.append ()
    b.append ()
    c.append ()
    d.append ()

# Define and solve the CVXPY problem.
x = cp.Variable(n)
# We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
soc_constraints = [cp.SOC (c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range (nn_h+nn_y)]
prob = cp.Problem(cp.Minimize(f.T@x),
                  soc_constraints + [F@x == g])
prob.solve()


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

