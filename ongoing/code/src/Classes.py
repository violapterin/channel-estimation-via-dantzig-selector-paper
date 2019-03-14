import numpy as np
import Constants as cst

class Ddss:
    def __init__(self, pp, y, gG):
        self.pp = pp
        self.y = y
        self.gG = gG

    def run(self):
        g = cp.Variable(n)
        m = cp.Variable(n)
        i = np.ones()
        socp_constr = (cvx.norm(A*x + b, 2) <= c.T*x + d)
        objective = cp.Minimize (i*m)
        constraints = [0 <= x, x <= 1]
        prob = cp.Problem(objective, constraints)
        self.sol = prob.solve()

class Oommpp:
    def __init__(self):
        pass

class Hh:
    def __init__(self):
        

class Ff_b:
    def __init__(self):
        np.random.normal();

class Ff_r:
    def __init__(self, len1, len2):
        
class Ww_b:
    def __init__(self, len1, len2):

class Ww_r:
    def __init__(self, len1, len2):
        
