

# import numpy as np
# import time

# a11, a12, a21, a22 = 1, 2, 3, 4
# b11, b12, b21, b22 = 5, 6, 7, 8


# A = np.array([[a11, a12], [a21, a22]])
# B = np.array([[b11, b12], [b21, b22]])


# start = time.time()
# for i in range(10000):
#     ans = A + B
# print("1は", time.time() - start)

# start = time.time()
# for i in range(10000):
#     ans = np.array([[a11+b11, a12+b12], [a21+b21, a22+b22]])
# print("2は", time.time() - start)



# # 関数の引数テスト

# def hoge(a, b, c, d):
#     print(a, b, c, d)
#     return a+b+c+d

# param = [
#     [1, 2, 3, 4],
#     [5, 6, 7, 8],
# ]

# hoge_list = [hoge(*param_list) for param_list in param]
# print(hoge_list)



# クラスのテスト

class HOGE:
    def __init__(self):
        self.a = 1
        self.b = 2

hoge = HOGE()
print(hoge.a)