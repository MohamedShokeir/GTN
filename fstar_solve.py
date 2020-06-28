#!/usr/bin/env python

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# material parameters :
f0 = 0.0035
fc = 0.04
fr = 0.056
q1 = 2
q2 = 1
delta = ((1./q1) - fc) / (fr - fc)

# polynominal parameters :
interv = 1.e-6
x1 = fc - interv
x2 = fc + interv
k1 = 1.
k2 = 0.
k3 = delta
k4 = fc * (1 - delta)

# A and b matrices :
A_ = np.array([[x1**5, x1**4, x1**3, x1**2, x1, 1.],
              [5.*x1**4, 4.*x1**3, 3.*x1**2, 2.*x1, 1., 0.],
              [20.*x1**3, 12.*x1**2, 6.*x1, 2., 0., 0.],
              [x2**5, x2**4, x2**3, x2**2, x2, 1.],
              [5.*x2**4, 4.*x2**3, 3.*x2**2, 2.*x2, 1., 0.],
              [20.*x2**3, 12.*x2**2, 6.*x2, 2., 0., 0.]])
print(np.linalg.cond(A_))

b_ = np.array([[k1*x1 + k2],
             [k1],
             [0.],
             [k3*x2 + k4],
             [k3],
             [0.]])

def main():

    def solve_coeff(A, b):
        """  
        Function that reutrns x, given A and b, in Ax = b
        """
        coeff = np.linalg.inv(A).dot(b)
        #coeff = linalg.solve(A, b)
        coeffs = coeff.flatten()
        np.savetxt("coeff_fstar_solve.res", coeff, delimiter=" ")
        return coeffs

    def fstar_poly(f):
        """
        Function that calculates fstar while : 
        1) f < x1 
        2) x1 <= f <= x2 (fifth-degree polynimal)  
        3) f > x2
        """
        if 0 <= f < x1:
            #print("i am f(x)")
            return k1*f + k2
        elif x1 <= f <= x2:
            #print("i am h(x)")
            return f * (f * (f * (f * (a * f + b) + c) + d) + e) + g
            #return a*f**5 + b*f**4 + c*f**3 + d*f**2 + e*f + g
        else:
            #print("i am g(x)")
            return k3*f + k4
   
    def dfstar_poly(f):
        """
        Function that calculates the derivative of fstar while : 
        1) f < x1 
        2) x1 <= f <= x2 (fifth-degree polynimal)  
        3) f > x2
        """
        if 0 <= f < x1:
            #print("i am f(x)")
            return k1 
        elif x1 <= f <= x2:
            #print("i am h(x)")
            return  f * (f * (f * (a * f * 5. + b * 4.) + c * 3.) + d * 2.) + e
            #return 5*a*f**4 + 4*b*f**3 + 3*c*f**2 + 2*d*f + e
        else:
            #print("i am g(x)")
            return k3
    
    def fstar_normal(f):
        """
        Function that calculates fstar as in the literature : 
        1) f* = f {while f < fc} 
        2) f* = fc + delta(f - fc) {otherwise}
        """
        if 0. <= f < fc:
            return k1*f + k2
        else:
            return k3*f + k4 

    # plot :
    a, b, c, d, e, g = solve_coeff(A_, b_)
    print(f"coeffs : {a, b, c, d, e, g}")
    #print(f"verif : {np.allclose(np.dot(A_, np.array([[a], [b], [c], [d], [e], [g]])), b_)}")
    #print(f"verif : {np.dot(A_, np.linalg.solve(A_, b_)) == b_}")
    f = np.linspace(0, 0.07, 1000)
    y_poly = [fstar_poly(i) for i in f]
    y_dpoly = [dfstar_poly(i) for i in f]
    y_normal = [fstar_normal(i) for i in f]
    # print(f"y_poly : {y_poly}") 
    # print(f"y_normal : {y_normal}") 
    
    plt.figure(figsize=(10, 8))
    plt.xlim(0, fr)
    plt.ylim(0, q2/2.)
    plt.plot(f, y_poly, label=r"$f_*$ polynominal")
    plt.plot(f, y_normal, label=r"$f_*$")
    plt.legend(loc='best')
    plt.savefig('fstar_solve.pdf', bbox_inches='tight')
    plt.show()
   
    plt.figure(figsize=(10, 8))
    plt.xlim(0, 0.06)
    #plt.ylim(0, 0.5)
    plt.plot(f, y_dpoly, '-o', label=r"$df_*$ polynominal")
    plt.legend(loc='best')
    plt.savefig('dfstar_solve.pdf', bbox_inches='tight')
    plt.show()


main()
