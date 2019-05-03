import matplotlib.pyplot as plt # plotting functions
import numpy as np
import constants as cst
#import os # getcwd
#import cvxpy as cp

import helpers as hp
import matrices as mt

#print (np.random.randint(7, size=10))
#print (np.sin (1J *np.pi)/ np.pi)

a =np.array ([1,2,4,5,7,8])
b1 =np.array ([1,-1,0,1,3,6])
b2 =np.array ([2,3,4,4.5,1,1])
b3 =np.array ([1,3,5,7,5,3])
#d_b1 =hp.Data(b1,"1st")
#d_b2 =hp.Data(b1,"2nd")
#d_b3 =hp.Data(b1,"3rd")
list_dat =[hp.Data(b1,"1st"), hp.Data(b2,"2nd"), hp.Data(b3,"3rd")]

hp.draw (a, list_dat, "Number", "Value", "Test")

#print (np.linalg.norm (b, ord='fro'))
#print(np.kron(a,b))
#print (type (np.reshape (a, (1,-1))[0]))
#print(np.reshape(np.kron(a,b), (1,-1)))

#aa=np.array ([[1,2,3], [4,5,6]])
#v=np.array ([1,2,3])
#print (np.shape(aa)[0])
#print (np.shape(v)[0])

#a =np.array ([[1,2],[4,5]])
#b =np.array ([1,-1])
#a[1][1]=3
#print (a @ b)
#def fun(v):
#    return 2*v
#v =np.array([1,2,3])
#print (fun(v))

#c=np.array([1,2,3,4,5,6])
#c=np.zeros(6)
#c=[1.5*i for i in range(6)]
#d=np.array([4,5,6,7,8,9])
#print (np.add (c, 1J *d))
#c_re =np.reshape(c, (2,-1))
#print (c_re[1])
#print (type(c_re[1]))
#print (np.array (np.reshape(c, (2,-1))[1]))
#print (c +1J *d)

#
# # # # # # # # # # # # # # # # # # # #
#
#list_x =[1, 2, 4]
#list_list_y =[[3.5, 6.2, 7], [2.1, 2.2, 2.3], [-0.2, 1, 3]]
#
#len_line =len(list_list_y)
#arr_x =np.array (list_x)
#
#for i in range (len_line):
#    arr_y_i =np.array(list_list_y[i])
#    plt.plot (arr_x, arr_y_i)
#
#out_fig_path =os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)) +"/plt/test.png"
##print (out_fig_path)
#plt.savefig (out_fig_path, bbox_inches ="tight")


