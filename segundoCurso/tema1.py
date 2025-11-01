import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def Dicotomia (a, b, fExpr, maxIter = 100, tol = 1.e-5, iter = [-1]):

    x = sp.symbols('x', real=True) 
    f = sp.Lambda(x,fExpr)

    xAprox = np.zeros(maxIter)

    if f(a) * f(b) > 0:
        print("La funcion no satisface las condiciones de Bolzano en el intervalo [a, b]")
        return None

    for k in range(0,maxIter):
        xAprox[k] = (a+b) / 2

        if f(xAprox[k]) == 0: break

        if f(a) * f(xAprox[k]) < 0:
            b = xAprox[k]
        else:
           a = xAprox[k]

        for value in iter:
            if value == k + 1:
                print('Iteración:', k+1, '| Aproximación:', round(xAprox[k], 10), '| Error:', '{:.10f}'.format(abs(a-b)/pow(2,k+1)), '| Error relativo:', '{:.10f}'.format(relativeError) )

        relativeError = np.abs( xAprox[k]-xAprox[k-1] ) / np.abs( xAprox[k] )
        if ( (k > 0) and ( relativeError < tol ) ): break

    print('Número de iteraciones realizadas: ', k+1) 
          # NOTA: Contamos 1 más, k+1, porque empezamos el bucle en 0
    print('Aproximación de la raíz: ', xAprox[k])

def InterpolacionLagrange(xValues, yValues):

    if len(xValues) != len(yValues):
        print("Los arrays de datos deben tener la misma longitud")
        return None

    x = sp.Symbol("x", real = True)

    xCoef = np.array([-2, 0, 1, 3])
    yCoef = np.array([-16, 0, -1, 9])

    n = len(xCoef)

    # Almacenaremos en "tabla" la matriz de diferencias divididas
    table = np.zeros([n, n])

    # La primera columna serán los datos en y
    table[:,0] = yCoef

    # Necesitamos un doble bucle para crear el resto de "tabla"
    for j in range(1,n):
        for i in range(n-j):
            table[i,j] = (table[i+1,j-1] - table[i,j-1]) / (xCoef[i+j]-xCoef[i])

    # Definimos la expresión para el Polinomio de Lagrange (versión Newton)
    pExpr = table[0,0]
    multiplica = sp.S('1')
    for k in range(1,n):
        multiplica = multiplica * (x - xCoef[k-1])
        pExpr = pExpr + table[0,k] * multiplica

    # Creamos la función lambdify para dibujarla
    return sp.lambdify(x,pExpr)

def calculateDomain(expr):
    x = sp.symbols('x', real=True)
    domain = sp.S.Reals

    # Get denominator of expr
    denom = expr.as_numer_denom()[1]

    # Exclude points where denominator is zero
    undefinedPoints = sp.solveset(denom, x, domain=sp.S.Reals)
    if undefinedPoints:
        domain = domain - undefinedPoints

    # Check for square roots of negative numbers
    for pow in expr.atoms(sp.Pow):
        base, exp = pow.as_base_exp()
        if exp == sp.Rational(1, 2):
            # Check if sqrt is in denominator (is a subexpression of denom)
            if denom.has(pow):
                # radicand must be strictly greater than 0 (denominator)
                domain = domain & sp.solveset(base > 0, x, domain=sp.S.Reals)
            else:
                # radicand must be greater or equal to 0
                domain = domain & sp.solveset(base >= 0, x, domain=sp.S.Reals)

    # Check for logarithms of non-positive numbers
    if expr.has(sp.log):
        logArgs = [arg.args[0] for arg in expr.atoms(sp.log)]
        for arg in logArgs:
            domain = domain & sp.solveset(arg > 0, x, domain=sp.S.Reals)

    return domain