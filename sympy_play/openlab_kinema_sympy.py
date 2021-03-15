"""openlabのロボットのいろいろ"""

import sympy as sy
from sympy import pi

# # qを定数とするとき
# q1, q2, q3, q4, q5, q6 = sy.symbols("q1, q2, q3, q4, q5, q6")
# qvec = sy.Matrix([[q1, q2, q3, q4, q5, q6]]).T

# l0, l1, l2, l3, l4, l5, l6, h45 = sy.symbols("l0, l1, l2, l3, l4, l5, l6, h45")

# posi_1 = sy.Matrix([[0, 0, 0]]).T
# posi_2 = posi_1 + sy.Matrix([[0, 0, l0]]).T
# posi_3 = posi_2 + sy.Matrix([[0, 0, l1]]).T
# posi_4 = posi_3 + sy.Matrix([[l2*sy.sin(q2)*sy.cos(q1)],
#                             [l2*sy.sin(q2)*sy.sin(q1)],
#                             [l2*sy.cos(q2)]])
# posi_5 = posi_4 + sy.Matrix([[l3*sy.sin(q2+q3)*sy.cos(q1)],
#                             [l3*sy.sin(q2+q3)*sy.sin(q1)],
#                             [l3*sy.cos(q2+q3)]])
# posi_6 = posi_5 + sy.Matrix([[(l4**2+h45**2)**(1/2)*sy.sin(q2+q3+q4+sy.atan(h45/l4))*sy.cos(q1)],
#                             [(l4**2+h45**2)**(1/2)*sy.sin(q2+q3+q4+sy.atan(h45/l4))*sy.sin(q1)],
#                             [(l4**2+h45**2)**(1/2)*sy.cos(q2+q3+q4+sy.atan(h45/l4))]])
# posi_7 = posi_6 + sy.Matrix([[l5*sy.sin(q2+q3+q4)*sy.cos(q1)],
#                             [l5*sy.sin(q2+q3+q4)*sy.sin(q1)],
#                             [l5*sy.cos(q2+q3+q4)]])

# A = q5
# B = pi/2 + q2 + q3 + q4
# C = q1
# D = q6

# def rotate_x(x):
#     return sy.Matrix([[1, 0, 0],
#                      [0, sy.cos(x), -sy.sin(x)],
#                      [0, sy.sin(x), sy.cos(x)]])

# def rotate_y(x):
#     return sy.Matrix([[sy.cos(x), 0, sy.sin(x)],
#                      [0, 1, 0],
#                      [-sy.sin(x), 0, sy.cos(x)]])

# def rotate_z(x):
#     return sy.Matrix([[sy.cos(x), -sy.sin(x), 0],
#                      [sy.sin(x), sy.cos(x), 0],
#                      [0, 0, 1]])

# posi_8 = posi_7 + rotate_z(C) * rotate_y(B) * rotate_x(A) * sy.Matrix([[l6*sy.cos(D), 0, l6*sy.sin(D)]]).T




# qを時間の関数とするとき
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
posi_4 = posi_3 + sy.Matrix([[l2*sy.sin(q2(t))*sy.cos(q1(t))],
                            [l2*sy.sin(q2(t))*sy.sin(q1(t))],
                            [l2*sy.cos(q2(t))]])
posi_5 = posi_4 + sy.Matrix([[l3*sy.sin(q2(t)+q3(t))*sy.cos(q1(t))],
                            [l3*sy.sin(q2(t)+q3(t))*sy.sin(q1(t))],
                            [l3*sy.cos(q2(t)+q3(t))]])
posi_6 = posi_5 + sy.Matrix([[(l4**2+h45**2)**(1/2)*sy.sin(q2(t)+q3(t)+q4(t)+sy.atan(h45/l4))*sy.cos(q1(t))],
                            [(l4**2+h45**2)**(1/2)*sy.sin(q2(t)+q3(t)+q4(t)+sy.atan(h45/l4))*sy.sin(q1(t))],
                            [(l4**2+h45**2)**(1/2)*sy.cos(q2(t)+q3(t)+q4(t)+sy.atan(h45/l4))]])
posi_7 = posi_6 + sy.Matrix([[l5*sy.sin(q2(t)+q3(t)+q4(t))*sy.cos(q1(t))],
                            [l5*sy.sin(q2(t)+q3(t)+q4(t))*sy.sin(q1(t))],
                            [l5*sy.cos(q2(t)+q3(t)+q4(t))]])

A = q5(t)
B = pi/2 + q2(t) + q3(t) + q4(t)
C = q1(t)
D = q6(t)

def rotate_x(x):
    return sy.Matrix([[1, 0, 0],
                     [0, sy.cos(x), -sy.sin(x)],
                     [0, sy.sin(x), sy.cos(x)]])

def rotate_y(x):
    return sy.Matrix([[sy.cos(x), 0, sy.sin(x)],
                     [0, 1, 0],
                     [-sy.sin(x), 0, sy.cos(x)]])

def rotate_z(x):
    return sy.Matrix([[sy.cos(x), -sy.sin(x), 0],
                     [sy.sin(x), sy.cos(x), 0],
                     [0, 0, 1]])

posi_8 = posi_7 + rotate_z(C) * rotate_y(B) * rotate_x(A) * sy.Matrix([[l6*sy.cos(D), 0, l6*sy.sin(D)]]).T




# print("posi_1 = ", sy.simplify(posi_1))
# print("posi_2 = ", sy.simplify(posi_2))
# print("posi_3 = ", sy.simplify(posi_3))
# print("posi_4 = ", sy.simplify(posi_4))
# print("posi_5 = ", sy.simplify(posi_5))
# print("posi_6 = ", sy.simplify(posi_6))
# print("posi_7 = ", sy.simplify(posi_7))
# print("posi_8 = ", sy.simplify(posi_8))



j1 = posi_1.jacobian(qvec)
j2 = posi_2.jacobian(qvec)
j3 = posi_3.jacobian(qvec)
j4 = posi_4.jacobian(qvec)
j5 = posi_5.jacobian(qvec)
j6 = posi_6.jacobian(qvec)
j7 = posi_7.jacobian(qvec)
j8 = posi_8.jacobian(qvec)

# print("j1 = ", sy.simplify(j1))
# print("j2 = ", sy.simplify(j2))
# print("j3 = ", sy.simplify(j3))
# print("j4 = ", sy.simplify(j4))
# print("j5 = ", sy.simplify(j5))
# print("j6 = ", sy.simplify(j6))
# print("j7 = ", sy.simplify(j7))
# print("j8 = ", sy.simplify(j8))

dj1 = sy.diff(j1, t)
dj2 = sy.diff(j2, t)
dj3 = sy.diff(j3, t)
dj4 = sy.diff(j4, t)
dj5 = sy.diff(j5, t)
dj6 = sy.diff(j6, t)
dj7 = sy.diff(j7, t)
dj8 = sy.diff(j8, t)

dj1 = sy.simplify(dj1)
dj2 = sy.simplify(dj2)
dj3 = sy.simplify(dj3)
dj4 = sy.simplify(dj4)
dj5 = sy.simplify(dj5)
dj6 = sy.simplify(dj6)
dj7 = sy.simplify(dj7)
dj8 = sy.simplify(dj8)

# dq1 = sy.Function("dq1")
# dq2 = sy.Function("dq2")
# dq3 = sy.Function("dq3")
# dq4 = sy.Function("dq4")
# dq5 = sy.Function("dq5")

dq1, dq2, dq3, dq4, dq5, dq6 = sy.symbols("dq1, dq2, dq3, dq4, dq5, dq6") 

dj1 = dj1.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
                (sy.Derivative(q6(t), t), dq6)])
dj2 = dj2.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
                (sy.Derivative(q6(t), t), dq6)])
dj3 = dj3.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
                (sy.Derivative(q6(t), t), dq6)])
dj4 = dj4.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
                (sy.Derivative(q6(t), t), dq6)])
dj5 = dj5.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
                (sy.Derivative(q6(t), t), dq6)])
dj6 = dj6.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
               (sy.Derivative(q6(t), t), dq6)])
dj7 = dj7.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
               (sy.Derivative(q6(t), t), dq6)])
dj8 = dj8.subs([(sy.Derivative(q1(t), t), dq1),
                (sy.Derivative(q2(t), t), dq2),
                (sy.Derivative(q3(t), t), dq3),
                (sy.Derivative(q4(t), t), dq4),
                (sy.Derivative(q5(t), t), dq5),
               (sy.Derivative(q6(t), t), dq6)])


# dj1 = dj1.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj2 = dj2.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj3 = dj3.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj4 = dj4.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj5 = dj5.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj6 = dj6.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])
# dj7 = dj7.subs([(q1(t), q1),
#                 (q2(t), q2),
#                 (q3(t), q3),
#                 (q4(t), q4),
#                 (q5(t), q5)])


#print("dj1 = ", dj1)
# print("dj2 = ", dj2)
# print("dj3 = ", dj3)
# print("dj4 = ", dj4)
# print("dj5 = ", dj5)
# print("dj6 = ", dj6)
# print("dj7 = ", dj7)
print("dj8 = ", dj8)

