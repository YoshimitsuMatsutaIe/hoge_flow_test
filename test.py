import baxter_utils_3
import baxter_utils_new

import numpy as np
import time


#q = np.zeros((7, 1))
q1_init = 0
q2_init = -31
q3_init = 0
q4_init = 43
q5_init = 0
q6_init = 72
q7_init = 0
q = np.array([[q1_init, q2_init, q3_init, q4_init, q5_init, q6_init, q7_init]]).T * np.pi / 180

start = time.time()
BaxterKinema = baxter_utils_3.BaxterKinematics3(
    L = 278e-3,
    h = 64e-3,
    H = 1104e-3,
    L0 = 270.35e-3,
    L1 = 69e-3,
    L2 = 364.35e-3,
    L3 = 69e-3,
    L4 = 374.29e-3,
    L5 = 10e-3,
    L6 = 368.3e-3
)
# for i in range(10000):
#     x = BaxterKinema.origins(q)
#     J = BaxterKinema.jacobi_all(q)
# print("3 = ", time.time() - start)

x = BaxterKinema.origins(q)
J = BaxterKinema.jacobi_all(q)
#print(x)
print(J, '\n')

start_2 = time.time()
BaxterKinema_2 = baxter_utils_new.BaxterKinematicsNew(
    L = 278e-3,
    h = 64e-3,
    H = 1104e-3,
    L0 = 270.35e-3,
    L1 = 69e-3,
    L2 = 364.35e-3,
    L3 = 69e-3,
    L4 = 374.29e-3,
    L5 = 10e-3,
    L6 = 368.3e-3
)
# for i in range(10000):
#     x_2 = BaxterKinema_2.origins(q)
#     J_2 = BaxterKinema_2.calc_jacobi(q)
# print("new = ", time.time() - start_2)


x_2 = BaxterKinema_2.origins(q)
J_2 = BaxterKinema_2.calc_jacobi(q)
#print(x_2)
print(J_2)