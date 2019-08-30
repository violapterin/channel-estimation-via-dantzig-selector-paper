#! /usr/bin/env python3

import matplotlib.pyplot as plt # plotting functions
import os # getcwd
import numpy as np
import cvxpy as cp

#import functions as fct
#import constants as cst

def Foo (a,b):
    print (a)
    print (b)

def Baz (fun, x):
    print ("call!")
    fun(x)

Baz (lambda x: Foo(x,2), 3)

quit()

for _ in range (500):
    print ("experiment:", _)
    pp =np.random.normal (0, 1, (10,50))
    y =np.random.normal (0, 1, (10))
    g = cp.Variable (50)
    prob = cp.Problem (
        cp.Minimize (cp.norm (g, 1)),
        [cp.norm (pp.conj().T @ (pp @ g - y), "inf") <= 0.5])
    prob.solve (solver = cp.CVXOPT)
    g_hat = g.value
    print (np.linalg.norm (g_hat, ord=1))





