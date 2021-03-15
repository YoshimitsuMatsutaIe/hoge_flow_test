"""baxterロボの同時変換行列などをsympyで式展開し整理
・左腕のみ
"""

import sympy as sy
import math

t = sy.Symbol("t")
# theta1, theta2, theta3, theta4, theta5, theta6, theta7 = sy.symbols("theta1, theta2, theta3, theta4, theta5, theta6, theta7")
# c1 = sy.cos(theta1)
# s1 = sy.sin(theta1)
# c2 = sy.cos(theta2)
# s2 = sy.sin(theta2)
# c3 = sy.cos(theta3)
# s3 = sy.sin(theta3)
# c4 = sy.cos(theta4)
# s4 = sy.sin(theta4)
# c5 = sy.cos(theta5)
# s5 = sy.sin(theta5)
# c6 = sy.cos(theta6)
# s6 = sy.sin(theta6)
# c7 = sy.cos(theta7)
# s7 = sy.sin(theta7)



# ### thetaを時間の関数にしたいとき ###
theta1 = sy.Function("theta1")
theta2 = sy.Function("theta2")
theta3 = sy.Function("theta3")
theta4 = sy.Function("theta4")
theta5 = sy.Function("theta5")
theta6 = sy.Function("theta6")
theta7 = sy.Function("theta7")
c1 = sy.cos(theta1(t))
s1 = sy.sin(theta1(t))
c2 = sy.cos(theta2(t))
s2 = sy.sin(theta2(t))
c3 = sy.cos(theta3(t))
s3 = sy.sin(theta3(t))
c4 = sy.cos(theta4(t))
s4 = sy.sin(theta4(t))
c5 = sy.cos(theta5(t))
s5 = sy.sin(theta5(t))
c6 = sy.cos(theta6(t))
s6 = sy.sin(theta6(t))
c7 = sy.cos(theta7(t))
s7 = sy.sin(theta7(t))


L, h, H = sy.symbols("L, h, H")
L0, L1, L2, L3, L4, L5, L6 = sy.symbols("L0, L1, L2, L3, L4, L5, L6")


# 同時変換行列

T_Wo_BL = sy.Matrix([[math.sqrt(2) / 2, math.sqrt(2) / 2, 0, L],
                     [-math.sqrt(2) / 2, math.sqrt(2) / 2, 0, -h],
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


# ジョイント角度ベクトル
q = sy.Matrix([[theta1(t), theta2(t), theta3(t), theta4(t), theta5(t), theta6(t), theta7(t)]]).T


jacobi_r_Wo_BL = r_Wo_BL.jacobian(q)
jacobi_r_Wo_0 = r_Wo_0.jacobian(q)
jacobi_r_Wo_1 = r_Wo_1.jacobian(q)
jacobi_r_Wo_2 = r_Wo_2.jacobian(q)
jacobi_r_Wo_3 = r_Wo_3.jacobian(q)
jacobi_r_Wo_4 = r_Wo_4.jacobian(q)
jacobi_r_Wo_5 = r_Wo_5.jacobian(q)
jacobi_r_Wo_6 = r_Wo_6.jacobian(q)
jacobi_r_Wo_7 = r_Wo_7.jacobian(q)
jacobi_r_Wo_GL = r_Wo_GL.jacobian(q)


#print("BL = ", jacobi_r_Wo_BL)
# print("0 = ", jacobi_r_Wo_0)
# print("1 = ", jacobi_r_Wo_1)
# print("2 = ", jacobi_r_Wo_2)
#print("3 = ", jacobi_r_Wo_3)
#print("4 = ", jacobi_r_Wo_4)
#print("5 = ", jacobi_r_Wo_5)
#print("6 = ", jacobi_r_Wo_6)
#print("7 = ", jacobi_r_Wo_7)
# #print("GL = ", jacobi_r_Wo_GL)

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


dtheta1 = sy.Function("dtheta1")
dtheta2 = sy.Function("dtheta2")
dtheta3 = sy.Function("dtheta3")
dtheta4 = sy.Function("dtheta4")
dtheta5 = sy.Function("dtheta5")
dtheta6 = sy.Function("dtheta6")
dtheta7 = sy.Function("dtheta7")


djacobi_BL = djacobi_BL.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                                (sy.Derivative(theta2(t), t), dtheta2(t)),
                                (sy.Derivative(theta3(t), t), dtheta3(t)),
                                (sy.Derivative(theta4(t), t), dtheta4(t)),
                                (sy.Derivative(theta5(t), t), dtheta5(t)),
                                (sy.Derivative(theta6(t), t), dtheta6(t)),
                                (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_0 = djacobi_0.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_1 = djacobi_1.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_2 = djacobi_2.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_3 = djacobi_3.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_4 = djacobi_4.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_5 = djacobi_5.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_6 = djacobi_6.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_7 = djacobi_7.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                            (sy.Derivative(theta2(t), t), dtheta2(t)),
                            (sy.Derivative(theta3(t), t), dtheta3(t)),
                            (sy.Derivative(theta4(t), t), dtheta4(t)),
                            (sy.Derivative(theta5(t), t), dtheta5(t)),
                            (sy.Derivative(theta6(t), t), dtheta6(t)),
                            (sy.Derivative(theta7(t), t), dtheta7(t)),])
djacobi_GL = djacobi_GL.subs([(sy.Derivative(theta1(t), t), dtheta1(t)),
                                (sy.Derivative(theta2(t), t), dtheta2(t)),
                                (sy.Derivative(theta3(t), t), dtheta3(t)),
                                (sy.Derivative(theta4(t), t), dtheta4(t)),
                                (sy.Derivative(theta5(t), t), dtheta5(t)),
                                (sy.Derivative(theta6(t), t), dtheta6(t)),
                                (sy.Derivative(theta7(t), t), dtheta7(t)),])

#print("BL = ", djacobi_BL)
# print("0 = ", djacobi_0)
# print("1 = ", djacobi_1)
#print("2 = ", djacobi_2)
#print("3 = ", djacobi_3)
#print("4 = ", djacobi_4)
#print("5 = ", sy.simplify(djacobi_5))
#print("6 = ", sy.simplify(djacobi_6))
#print("7 = ", sy.simplify(djacobi_7))
print("GL = ", sy.simplify(djacobi_GL))

