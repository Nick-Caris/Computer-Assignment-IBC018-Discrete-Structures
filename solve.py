from sympy import symbols

from toolbox import import_function
import sympy as sy

#
# function, initialStates = import_function.importFunction('FILENAME')
# Step 1
# function =
#
# Step 2
# characteristicEquation =
#
# Step 3
# FindRoot =
#
# Step 4
# General Solution =
#
# Step 5
# ExactValue =
#
# directFormula =
#
# toolbox_import.filename(directFormula)
# Job test

# x, y, z = symbols('x y z')
# expression1 = x * 4 ** 6 + y
# expression2 = x
# solution = sy.solve(sy.Eq(expression1, expression2), particular=True)
# print(solution)
# expression2 = expression2.subs(x, solution[0][x])
# expression1 = expression1.subs(x, solution[0][x])
# print(expression1)
# print(expression2)
# solution = sy.solve(sy.Eq(expression1, expression2), particular=True)
# print(solution)

A, B, n = symbols('A B n')
expression3 = A * n + B
expression4 = 2 * (A * (n - 1) + B) + n + 5
print(sy.solve(sy.Eq(expression3, expression4)))
print(sy.solve(sy.Eq(expression3, expression4), [A,B], particular=True))
print(range(3 + 1, 0))
for i in range(3 + 1, 0, -1):
    print(i)