import sympy as sy




q, dq, q_min, q_max, sigma, alpha_u, alpha_l = sy.symbols('q, dq, q_min, q_max, sigma, alpha_u, alpha_l')





# def jl_alpha_upper(dq, sigma):
#     return 1 - sy.exp(-(max(dq, 0) ** 2) / (2 * sigma ** 2))

# def jl_alpha_lower(dq, sigma):
#     return 1 - sy.exp(-(min(dq, 0) ** 2) / (2 * sigma ** 2))

def jl_s(q, q_min, q_max):
    return (q - q_min) / (q_max - q_min)

def jl_b(q, dq, q_min, q_max, sigma, alpha_u, alpha_l):
    """???を計算"""
    s = jl_s(q, q_min, q_max)
    d = 4 * s * (1 - s)
    return s * (alpha_u * d + (1 - alpha_u)) + (1 - s) * (alpha_l * d + (1 - alpha_l))


b = jl_b(q, dq, q_min, q_max, sigma, alpha_u, alpha_l)

a = b ** (-2)

da = sy.diff(a, q)

da = sy.simplify(da)

print(da)