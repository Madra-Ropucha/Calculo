import sympy as sp

x, y = sp.symbols('x:y', real=True)

p3 = sp.plot_implicit(sp.Or(y > x, y < -2 * x))

# x = sp.Symbol('x')
# y = sp.Symbol('y')

# f = 3 * x + 5 * y
# print(f)

# pi = sp.pi

# print(sp.sin(pi))
# print(sp.cos(pi))

# print(sp.log(sp.E))

# print(sp.exp(2))

# a,s = sp.symbols('x,y', real=True)

# expr = (x-3)*(x-3)**2*(y-2)

# print(expr.factor())

# expr = sp.exp(x+1)-5
# print(sp.solve(expr))

# q = sp.symbols('x', real=True)
# w = sp.plot(sp.sin(x), sp.cos(x), (x, 0, 4*sp.pi), show=False)
# w[0].line_color = 'r'
# w[1].line_color = 'b'
# w.xlabel = 'x'
# w.ylabel = 'y'

# w.legend = True
# w.show()