"""sympyでrmp(from GDS)を計算"""
import sympy as sy
import time


start = time.time()

x0, y0, z0 = sy.symbols('x0, y0, z0')
dx0, dy0, dz0 = sy.symbols('dx0, dy0, dz0')

x = sy.Matrix([[x0, y0, z0]]).T
dx = sy.Matrix([[dx0, dy0, dz0]]).T

#x_norm = x.norm(ord = 2)
x_norm = (x0**2 + y0**2 + z0**2) ** (1/2)
x_hat = x / x_norm

#print(x_norm)

sigma_alpha, sigma_gamma, w_u, w_l, alpha, epsilon = sy.symbols('sigma_alpha, sigma_gamma, w_u, w_l, alpha, epsilon')


alpha_x = sy.exp(-((x_norm**2) / (2 * sigma_alpha ** 2)))
gamma_x = sy.exp(-((x_norm**2) / (2 * sigma_gamma ** 2)))
w_x = gamma_x * w_u + (1 - gamma_x) * w_l
nabla_x_phi = (1 - sy.exp(-2 * alpha * x_norm)) / (1 + sy.exp(-2 * alpha * x_norm)) * x_hat


M = w_x  * ((1 - alpha_x) * nabla_x_phi * nabla_x_phi.T + (alpha_x + epsilon) * sy.eye(3))

print(sy.simplify(M))


# # 曲率校xi_Mの前半を計算
# partialx_m1 = M[:, 0:1].jacobian(x)
# partialx_m2 = M[:, 1:2].jacobian(x)
# partialx_m3 = M[:, 2:3].jacobian(x)

# # print(sy.simplify(partialx_m1))

# partialx_m1_dx = partialx_m1 * dx
# partialx_m2_dx = partialx_m2 * dx
# partialx_m3_dx = partialx_m3 * dx

# xi_M_before_before = partialx_m1_dx.row_join(partialx_m2_dx)
# xi_M_before_before = xi_M_before_before.row_join(partialx_m3_dx)

# xi_M_before = xi_M_before_before * dx

# #print(xi_M_before)

# # 曲率校xi_Mの後半を計算
# xi_M_after_before = dx.T * M * dx
# xi_M_after = (xi_M_after_before.jacobian(x)).T #　sympyにはナブラあるかも？

# # 合体
# xi_M = xi_M_before + xi_M_after
# xi_M = sy.simplify(xi_M)  # Elisじゃ無理
# print(xi_M)

print("計算時間 ", time.time() - start)