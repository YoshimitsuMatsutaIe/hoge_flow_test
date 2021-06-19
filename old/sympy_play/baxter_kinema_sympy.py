"""baxterロボの同時変換行列などをsympyで式展開し整理
・左腕のみ
"""

import sympy as sy
import math


# q1, q2, q3, q4, q5, q6, q7 = sy.symbols("q1, q2, q3, q4, q5, q6, q7")
# c0 = sy.cos(sy.pi / 4)
# s0 = sy.sin(sy.pi / 4)
# c1 = sy.cos(q1)
# s1 = sy.sin(q1)
# c2 = sy.cos(q2)
# s2 = sy.sin(q2)
# c3 = sy.cos(q3)
# s3 = sy.sin(q3)
# c4 = sy.cos(q4)
# s4 = sy.sin(q4)
# c5 = sy.cos(q5)
# s5 = sy.sin(q5)
# c6 = sy.cos(q6)
# s6 = sy.sin(q6)
# c7 = sy.cos(q7)
# s7 = sy.sin(q7)
# # ジョイント角度ベクトル
# q = sy.Matrix([[q1, q2, q3, q4, q5, q6, q7]]).T



### qを時間の関数にしたいとき ###
t = sy.Symbol("t")
c0 = sy.cos(sy.pi / 4)
s0 = sy.sin(sy.pi / 4)
q1 = sy.Function("q1")
q2 = sy.Function("q2")
q3 = sy.Function("q3")
q4 = sy.Function("q4")
q5 = sy.Function("q5")
q6 = sy.Function("q6")
q7 = sy.Function("q7")
c1 = sy.cos(q1(t))
s1 = sy.sin(q1(t))
c2 = sy.cos(q2(t))
s2 = sy.sin(q2(t))
c3 = sy.cos(q3(t))
s3 = sy.sin(q3(t))
c4 = sy.cos(q4(t))
s4 = sy.sin(q4(t))
c5 = sy.cos(q5(t))
s5 = sy.sin(q5(t))
c6 = sy.cos(q6(t))
s6 = sy.sin(q6(t))
c7 = sy.cos(q7(t))
s7 = sy.sin(q7(t))

# ジョイント角度ベクトル
q = sy.Matrix([[q1(t), q2(t), q3(t), q4(t), q5(t), q6(t), q7(t)]]).T

L, h, H = sy.symbols("L, h, H")
L0, L1, L2, L3, L4, L5, L6 = sy.symbols("L0, L1, L2, L3, L4, L5, L6")


# 同時変換行列

# 直書き
# T_Wo_BL = sy.Matrix([[math.sqrt(2) / 2, math.sqrt(2) / 2, 0, L],
#                      [-math.sqrt(2) / 2, math.sqrt(2) / 2, 0, -h],
#                      [0, 0, 1, H],
#                      [0, 0, 0, 1]])

# 直書きじゃない
T_Wo_BL = sy.Matrix([[c0, s0, 0, L],
                     [-s0, c0, 0, -h],
                     [0, 0, 1, H],
                     [0, 0, 0, 1]])

T_BL_0 = sy.Matrix([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, L0],
                    [0, 0, 0, 1]])

T_0_1 = sy.Matrix([[c1, -s1, 0, 0],
                   [s1, c1, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])

T_1_2 = sy.Matrix([[-s2, -c2, 0, L1],
                   [0, 0, 1, 0],
                   [-c2, s2, 0, 0],
                   [0, 0, 0, 1]])

T_2_3 = sy.Matrix([[c3, -s3, 0, 0],
                   [0, 0, -1, -L2],
                   [s3, c3, 0, 0],
                   [0, 0, 0, 1]])

T_3_4 = sy.Matrix([[c4, -s4, 0, L3],
                   [0, 0, 1, 0],
                   [-s4, -c4, 0, 0],
                   [0, 0, 0, 1]])

T_4_5 = sy.Matrix([[c5, -s5, 0, 0],
                   [0, 0, -1, -L4],
                   [s5, c5, 0, 0],
                   [0, 0, 0, 1]])

T_5_6 = sy.Matrix([[c6, -s6, 0, L5],
                   [0, 0, 1, 0],
                   [-s6, -c6, 0, 0],
                   [0, 0, 0, 1]])

T_6_7 = sy.Matrix([[c7, -s7, 0, 0],
                   [0, 0, -1, 0],
                   [s7, c7, 0, 0],
                   [0, 0, 0, 1]])

T_7_GL = sy.Matrix([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, L6],
                    [0, 0, 0, 1]])

T_0_7 = T_0_1 * T_1_2 * T_2_3 * T_3_4 * T_4_5 * T_5_6 * T_6_7
T_0_6 = T_0_1 * T_1_2 * T_2_3 * T_3_4 * T_4_5 * T_5_6
T_0_5 = T_0_1 * T_1_2 * T_2_3 * T_3_4 * T_4_5
T_0_4 = T_0_1 * T_1_2 * T_2_3 * T_3_4
T_0_3 = T_0_1 * T_1_2 * T_2_3
T_0_2 = T_0_1 * T_1_2


## 正しい
T_Wo_BL = T_Wo_BL
T_Wo_0 = T_Wo_BL * T_BL_0
T_Wo_1 = T_Wo_BL * T_BL_0 * T_0_1
T_Wo_2 = T_Wo_BL * T_BL_0 * T_0_2
T_Wo_3 = T_Wo_BL * T_BL_0 * T_0_3
T_Wo_4 = T_Wo_BL * T_BL_0 * T_0_4
T_Wo_5 = T_Wo_BL * T_BL_0 * T_0_5
T_Wo_6 = T_Wo_BL * T_BL_0 * T_0_6
T_Wo_7 = T_Wo_BL * T_BL_0 * T_0_7
T_Wo_GL = T_Wo_BL * T_BL_0 * T_0_7 * T_7_GL

# print("r_bar_Wo_BL = ", T_Wo_BL[0:3, 3:4])
# print("r_bar_Wo_0 = ", T_Wo_0[0:3, 3:4])
# print("r_bar_Wo_1 = ", T_Wo_1[0:3, 3:4])
# print("r_bar_Wo_2 = ", T_Wo_2[0:3, 3:4])
# print("r_bar_Wo_3 = ", T_Wo_3[0:3, 3:4])
# print("r_bar_Wo_4 = ", T_Wo_4[0:3, 3:4])
# print("r_bar_Wo_5 = ", T_Wo_5[0:3, 3:4])
# print("r_bar_Wo_6 = ", T_Wo_6[0:3, 3:4])
# print("r_bar_Wo_7 = ", T_Wo_7[0:3, 3:4])
# print("r_bar_Wo_GL = ", T_Wo_GL[0:3, 3:4])

#print("orig = ", T_Wo_GL[0:1, 3:4])
#print("simp = ", sy.simplify(T_Wo_GL[0:1, 3:4]))

# 世界座標系から見た局所座標系の原点位置
r_Wo_BL = T_Wo_BL[0:3, 3:4]
r_Wo_0 = T_Wo_0[0:3, 3:4]
r_Wo_1 = T_Wo_1[0:3, 3:4]
r_Wo_2 = T_Wo_2[0:3, 3:4]
r_Wo_3 = T_Wo_3[0:3, 3:4]
r_Wo_4 = T_Wo_4[0:3, 3:4]
r_Wo_5 = T_Wo_5[0:3, 3:4]
r_Wo_6 = T_Wo_6[0:3, 3:4]
r_Wo_7 = T_Wo_7[0:3, 3:4]
r_Wo_GL = T_Wo_GL[0:3, 3:4]

r_Wo_BL = sy.simplify(r_Wo_BL)
r_Wo_0 = sy.simplify(r_Wo_0)
r_Wo_1 = sy.simplify(r_Wo_1)
r_Wo_2 = sy.simplify(r_Wo_2)
r_Wo_3 = sy.simplify(r_Wo_3)
r_Wo_4 = sy.simplify(r_Wo_4)
r_Wo_5 = sy.simplify(r_Wo_5)
r_Wo_6 = sy.simplify(r_Wo_6)
r_Wo_7 = sy.simplify(r_Wo_7)
r_Wo_GL = sy.simplify(r_Wo_GL)

# print("BL = ", r_Wo_BL)
# print("0 = ", r_Wo_0)
# print("1 = ", r_Wo_1)
# print("2 = ", r_Wo_2)
# print("3 = ", r_Wo_3)
# print("4 = ", r_Wo_4)
# print("5 = ", r_Wo_5)
# print("6 = ", r_Wo_6)
# print("7 = ", r_Wo_7)
# print("GL = ", r_Wo_GL)


jacobi_r_Wo_BL = sy.simplify(r_Wo_BL.jacobian(q))
jacobi_r_Wo_0 = sy.simplify(r_Wo_0.jacobian(q))
jacobi_r_Wo_1 = sy.simplify(r_Wo_1.jacobian(q))
jacobi_r_Wo_2 = sy.simplify(r_Wo_2.jacobian(q))
jacobi_r_Wo_3 = sy.simplify(r_Wo_3.jacobian(q))
jacobi_r_Wo_4 = sy.simplify(r_Wo_4.jacobian(q))
jacobi_r_Wo_5 = sy.simplify(r_Wo_5.jacobian(q))
jacobi_r_Wo_6 = sy.simplify(r_Wo_6.jacobian(q))
jacobi_r_Wo_7 = sy.simplify(r_Wo_7.jacobian(q))
jacobi_r_Wo_GL = sy.simplify(r_Wo_GL.jacobian(q))

# print("BL = ", jacobi_r_Wo_BL)
# print("0 = ", jacobi_r_Wo_0)
# print("1 = ", jacobi_r_Wo_1)
# print("2 = ", jacobi_r_Wo_2)
# print("3 = ", jacobi_r_Wo_3)
# print("4 = ", jacobi_r_Wo_4)
# print("5 = ", jacobi_r_Wo_5)
# print("6 = ", jacobi_r_Wo_6)
# print("7 = ", jacobi_r_Wo_7)
# print("GL = ", jacobi_r_Wo_GL)



djacobi_BL = sy.diff(jacobi_r_Wo_BL, t)
djacobi_0 = sy.diff(jacobi_r_Wo_0, t)
djacobi_1 = sy.diff(jacobi_r_Wo_1, t)
djacobi_2 = sy.diff(jacobi_r_Wo_2, t)
djacobi_3 = sy.diff(jacobi_r_Wo_3, t)
djacobi_4 = sy.diff(jacobi_r_Wo_4, t)
djacobi_5 = sy.diff(jacobi_r_Wo_5, t)
djacobi_6 = sy.diff(jacobi_r_Wo_6, t)
djacobi_7 = sy.diff(jacobi_r_Wo_7, t)
djacobi_GL = sy.diff(jacobi_r_Wo_GL, t)


# # dqを時間の変数にしたいとき
# dq1 = sy.Function("dq1")
# dq2 = sy.Function("dq2")
# dq3 = sy.Function("dq3")
# dq4 = sy.Function("dq4")
# dq5 = sy.Function("dq5")
# dq6 = sy.Function("dq6")
# dq7 = sy.Function("dq7")

# djacobi_BL = djacobi_BL.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                                 (sy.Derivative(q2(t), t), dq2(t)),
#                                 (sy.Derivative(q3(t), t), dq3(t)),
#                                 (sy.Derivative(q4(t), t), dq4(t)),
#                                 (sy.Derivative(q5(t), t), dq5(t)),
#                                 (sy.Derivative(q6(t), t), dq6(t)),
#                                 (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_0 = djacobi_0.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_1 = djacobi_1.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_2 = djacobi_2.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_3 = djacobi_3.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_4 = djacobi_4.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_5 = djacobi_5.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_6 = djacobi_6.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_7 = djacobi_7.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                             (sy.Derivative(q2(t), t), dq2(t)),
#                             (sy.Derivative(q3(t), t), dq3(t)),
#                             (sy.Derivative(q4(t), t), dq4(t)),
#                             (sy.Derivative(q5(t), t), dq5(t)),
#                             (sy.Derivative(q6(t), t), dq6(t)),
#                             (sy.Derivative(q7(t), t), dq7(t)),])
# djacobi_GL = djacobi_GL.subs([(sy.Derivative(q1(t), t), dq1(t)),
#                                 (sy.Derivative(q2(t), t), dq2(t)),
#                                 (sy.Derivative(q3(t), t), dq3(t)),
#                                 (sy.Derivative(q4(t), t), dq4(t)),
#                                 (sy.Derivative(q5(t), t), dq5(t)),
#                                 (sy.Derivative(q6(t), t), dq6(t)),
#                                 (sy.Derivative(q7(t), t), dq7(t)),])


# dqを時間に依らないとしたいとき
dq1 ,dq2, dq3, dq4, dq5, dq6, dq7 = sy.symbols('dq1 ,dq2, dq3, dq4, dq5, dq6, dq7')

djacobi_BL = djacobi_BL.subs([(sy.Derivative(q1(t), t), dq1),
                                (sy.Derivative(q2(t), t), dq2),
                                (sy.Derivative(q3(t), t), dq3),
                                (sy.Derivative(q4(t), t), dq4),
                                (sy.Derivative(q5(t), t), dq5),
                                (sy.Derivative(q6(t), t), dq6),
                                (sy.Derivative(q7(t), t), dq7)])
djacobi_0 = djacobi_0.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_1 = djacobi_1.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_2 = djacobi_2.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_3 = djacobi_3.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_4 = djacobi_4.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_5 = djacobi_5.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_6 = djacobi_6.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_7 = djacobi_7.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])
djacobi_GL = djacobi_GL.subs([(sy.Derivative(q1(t), t), dq1),
                            (sy.Derivative(q2(t), t), dq2),
                            (sy.Derivative(q3(t), t), dq3),
                            (sy.Derivative(q4(t), t), dq4),
                            (sy.Derivative(q5(t), t), dq5),
                            (sy.Derivative(q6(t), t), dq6),
                            (sy.Derivative(q7(t), t), dq7),])



djacob_BL = sy.simplify(djacobi_BL)
djacob_0 = sy.simplify(djacobi_0)
djacob_1 = sy.simplify(djacobi_1)
djacob_2 = sy.simplify(djacobi_2)
djacob_3 = sy.simplify(djacobi_3)
djacob_4 = sy.simplify(djacobi_4)
djacob_5 = sy.simplify(djacobi_5)
djacob_6 = sy.simplify(djacobi_6)
djacob_7 = sy.simplify(djacobi_7)
djacob_GL = sy.simplify(djacobi_GL)


print("BL = ", djacobi_BL)
print("0 = ", djacobi_0)
print("1 = ", djacobi_1)
print("2 = ", djacobi_2)
print("3 = ", djacobi_3)
print("4 = ", djacobi_4)
print("5 = ", djacobi_5)
print("6 = ", djacobi_6)
print("7 = ", djacobi_7)
print("GL = ", djacobi_GL)

