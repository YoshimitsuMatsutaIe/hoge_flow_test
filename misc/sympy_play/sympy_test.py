import sympy as sy
from sympy.printing.pycode import pycode


t = sy.Symbol("t")
q1 = sy.Function("q1")
q2 = sy.Function("q2")
q3 = sy.Function("q3")
q4 = sy.Function("q4")
q5 = sy.Function("q5")
q6 = sy.Function("q6")
qvec = sy.Matrix([[q1(t), q2(t), q3(t), q4(t), q5(t), q6(t)]]).T


l0, l1, l2, l3, l4, l5, l6, h45 = sy.symbols("l0, l1, l2, l3, l4, l5, l6, h45")

posi_1 = sy.Matrix([[0, 0, 0]]).T
posi_2 = posi_1 + sy.Matrix([[0, 0, l0]]).T
posi_3 = posi_2 + sy.Matrix([[0, 0, l1]]).T


code = pycode(posi_3)
print(code)