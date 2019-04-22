import matplotlib.pyplot as plt # plotting functions
import numpy as np
import constants as cst
#import os # getcwd
#import cvxpy as cp

#a=[[1,2,3,4],[5,6,7,8]]
#b=[[1,-1],[0,1]]
#print(np.reshape(np.kron(a,b), (1,-1)))

a =[1.3*i for i in range (5)]
print (type (a))

c=np.array([1,2,3,4,5,6])
#c=np.zeros(6)
#c=[1.5*i for i in range(6)]
#d=np.array([4,5,6,7,8,9])
#print (np.add (c, 1J *d))
#c_re =np.reshape(c, (2,-1))
#print (c_re[1])
#print (type(c_re[1]))
#print (np.array (np.reshape(c, (2,-1))[1]))
#print (c +1J *d)

#y=[1+2j,2,3,4]
#z1=[(y[2*i]) for i in range(2)]
#z2=[(y[2*i+1]) for i in range(2)]
#print(z1)
#print(z2)

#A=np.array ([[1,2],[3,4],[5,6]])
#B=np.array ([[1,2],[4,5]])
#b=np.array ([1,2])
#c=np.array ([[1],[2]])
#d=b.T
#
#print(A@b)
#print(A@c)
#print(A@d)
#print(A@b.T)

#
# # # # # # # # # # # # # # # # # # # #
#
#lst_x =[1, 2, 4]
#lst_lst_y =[[3.5, 6.2, 7], [2.1, 2.2, 2.3], [-0.2, 1, 3]]
#
#len_line =len(lst_lst_y)
#arr_x =np.array (lst_x)
#
#for i in range (len_line):
#    arr_y_i =np.array(lst_lst_y[i])
#    plt.plot (arr_x, arr_y_i)
#
#out_fig_path =os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)) +"/plt/test.png"
##print (out_fig_path)
#plt.savefig (out_fig_path, bbox_inches ="tight")


