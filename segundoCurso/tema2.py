import numpy as np
import sympy as sp

def NewtonRaphson(fExpr, x0, maxIter=100, tol=1.e-9):
    x = sp.symbols('x', real=True) # define la variable simbólica x

    fDerExpr = sp.diff(fExpr,x)

    f = sp.lambdify(x,fExpr)
    fDer = sp.lambdify(x,fDerExpr)

    xAprox = np.zeros(maxIter)
    xAprox[0] = x0

    for k in range(1,maxIter):
        if ( np.abs( fDer(xAprox[k-1]) ) < 1.e-14 ): break

        xAprox[k] = xAprox[k-1] - f(xAprox[k-1])/fDer(xAprox[k-1])

        if ( (k > 0) and (np.abs(xAprox[k]-xAprox[k-1]) / np.abs(xAprox[k]) < tol) ): break

    print('Número de iteraciones realizadas: ', k) 
    print('Aproximación de la raíz: ', xAprox[k])
