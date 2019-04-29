import numpy as np
import constants as cst


class Hh:
    def __init__(self):
        val =np.zeros((self.cst.nn_h, self.cst.nn_h))
    
    def refresh(self):
        for l in range (self.cst.ll)
            pG_l =(d_array /cst.lG) * np.sin (np.random.uniform(0, 2 * np.pi))
            tG_l =(d_array /cst.lG) * np.sin (np.random.uniform(0, 2 * np.pi))
            (arr_response (pG_l)).H @ arr_response (pG_l)
        
        for i in range(cst.nn_h)
            for j in range(cst.nn_h)
                val[i][j] =

class Ff_b:
    def __init__(self):
        mat =np.zeros ((cst.nn_h, cst.nn_h))


        np.random.normal (cst.normal_mean, cst.normal_std)

class Ff_r:
    def __init__(self, cst):
        self.cst =cst
        mat =np.zeros((self.cst.nn_h, self.cst.nn_h))

    def generate(self):
        for i in range(cst.nn_r)
            for j in range(cst.nn_y)
                idx_phase =np.random.randint (cst.num_quantization)
                matrix[i][j] =np.exp (2 * np.pi * 1J * idx_phase /num_quantization)
        
class Ww_b:
    def __init__(self, len1, len2):

class Ww_r:
    def __init__(self, len1, len2):
        
