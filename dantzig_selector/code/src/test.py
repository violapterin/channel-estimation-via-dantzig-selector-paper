import matplotlib.pyplot as plt # plotting functions
import numpy as np
import constants as cst
import os # getcwd
#import cvxpy as cp

import functions as fct


tt =np.array (range(3))

pp =np. array ([[4,5,6,3], [0,-2,7,4]])
g =np. array ([4,4,5,5])
y =pp @ g
ss =set ([0,3])
r =y

t =np.argmin (abs (pp [:, tt].conj().T @ r))
ss.add(t)
print (np. array (list(ss))*2)

#i =[1,2,4,5]
#aa =np.array ([[5,8,-2,-3,0,9], [3,6,7,9,10,11]])
#bb =aa [:, i]
#print (bb)

#a =[3,1,2]
#print (a)
#a.sort()
#print (a)

#print (np.random.randint(7, size=10))
#print (np.sin (1J *np.pi)/ np.pi)

#zz =mt.Zz()
#zz.generate()
#hh =mt.Hh()
#hh.generate()

#a =np.array ([[1+2J, 3+4.4J], [-1J, -2.3+3J]])
#print (a.conj().T)

#a = np.array ([ [(i+j^2) for i in range(5)] for j in range(5)])
#a = [ [(i+j^2) for i in range(5)] for j in range(5)]
#print (a)
#s =set([0,2,4])
#b =a [:, list(s)]
#print (b)


#a =np.array ([1,2,4,5,7,8])
#b1 =np.array ([1,-1,0,1,3,6])
#b2 =np.array ([2,3,4,4.5,1,1])
#b3 =np.array ([1,3,5,7,5,3])
#list_dat =[hp.Data(b1,"1st"), hp.Data(b2,"2nd"), hp.Data(b3,"3rd")]
#hp.draw (a, list_dat, "Number", "Value", "Test")

#print (np.linalg.norm (b, ord='fro'))
#print(np.kron(a,b))
#print (type (np.reshape (a, (1,-1))[0]))
#print(np.reshape(np.kron(a,b), (1,-1)))


