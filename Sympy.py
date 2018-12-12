from sympy import symbols, init_printing
import sympy as sp

# 2*sqrt(x)+5^y+(2/z) represented as a sympy expression
x, y, z = symbols('x y z')
expression = 2 * sp.sqrt(x) + 5 ** y + (2 / z)
print(expression)

# Add 1 to the expression
expression = expression + 1
print(expression)

# Subtract (2/z) from the expression (should resolve (2/z) - (2/z) to 0 automatically)
expression = expression - (2 / z)
print(expression)

# This does not automatically resolve:
expression = x * expression
print(expression)

# This expands the formula
expanded_expression = sp.expand(expression)
print(expanded_expression)

# Solves an equation (implicitly adds = 0 to an expression
print(sp.solve(expression))
print(sp.solve(x ** 2 + 2 * x - 10))
