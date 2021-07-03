
# from autograd import grad as nabla
# import autograd
# import autograd.numpy as agnp
# import numpy as np

# # def f(x):
# #     Q = np.array([[3, 2],
# #         [2, 7]])
# #     #return 1./2 * np.dot(x, np.dot(Q, x))
# #     return 1./2 * x @ Q @ x

# # x = np.array([1.0, 2.0])
# # print (nabla(f)(x))


# def M_stretch(x):
#     x = x.reshape(3, 1)
#     norm_x = agnp.linalg.norm(x)
#     hat_x = x / norm_x
#     alpha = agnp.exp(-norm_x**2 / (2*1**2))
#     gamma = agnp.exp(-norm_x**2 / (2*1**2))
#     w = gamma * 1 + (1 - gamma) * 1
    
#     def s_alpha(alpha, norm_x):
#         return (1 - agnp.exp(-2*alpha*norm_x)) / (1 + agnp.exp(-2*alpha*norm_x))
    
#     nabla_x_potential = s_alpha(1, norm_x) * hat_x
#     M = w * ((1 - alpha) * (nabla_x_potential @ nabla_x_potential.T) + \
#         (alpha + 1) * np.eye(3))
#     return M

# def partial_mi_s(x):
    
#     #print(M_stretch(x))
#     #x_ = np([[x[0,0], x[1,0], x[2, 0]]])
#     def m1(x):
        
#         return M_stretch(x)[:, 0]
    
#     def m2(x):
#         return M_stretch(x)[:, 1]
    
#     def m3(x):
#         return M_stretch(x)[:, 2]
    
#     partial_x_m1 = autograd.jacobian(m1)(x)
#     partial_x_m2 = autograd.jacobian(m2)(x)
#     partial_x_m3 = autograd.jacobian(m3)(x)
    
#     return [partial_x_m1, partial_x_m2, partial_x_m3]



# x = np.array([1., 1., 1.]).T
# dx = np.array([[1.0, 1.0, 1.0]]).T

# partial_ms = partial_mi_s(x)

# [print(j) for j in partial_ms]

# xi_1_ = [par @ dx for par in partial_ms]
# xi_1 = np.concatenate(xi_1_, axis = 1) @ dx

import numpy as np

a = np.array([[1, 2, 3]])
b, c, d = a
print(b, c, d)