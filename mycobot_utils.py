"""mycobotの色々な関数とクラス
・順運動学写像，ヤコビ行列など
"""

import numpy as np
from math import cos, sin, tan, pi, sqrt


class MyCobotKinematics:
    """ローカル座標系の原点，ヤコビ行列等を計算
    ・ちょっとはやい
    """
    
    def __init__(self, **kwargs):
        self.d1 = kwargs.pop('d1')
        self.a2 = kwargs.pop('a2')
        self.a3 = kwargs.pop('a3')
        self.d4 = kwargs.pop('d4')
        self.d5 = kwargs.pop('d5')
        self.d6 = kwargs.pop('d6')
    
    def homo_transf_matrix(self, theta, a, alpha, d):
        """DH記法で使う
        ・未完成
        """
        return np.array([
            [cos(theta), -sin(theta), 0, a],
            [sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha)],
            [sin(theta)*sin(alpha), cos(theta)*sin(alpha), -cos(alpha), d*cos(alpha)],
            [0, 0, 0, 1],
        ])
    
    
    def init_homo_transf_mat(self, q):
        """同時変換行列を作成
        ・
        """
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        DHparam = [
            [q1, 0, pi/2, self.d1],
            [q2, -self.a2, 0, 0],
            [q3, -self.a3, 0, 0],
            [q4, 0, pi/2, self.d4],
            [q5, 0, -pi/2, self.d5],
            [q6, 0, 0, self.d6],
        ]
        
        # 個別
        self.T_i_j = [self.homo_transf_matrix(*param) for param in DHparam]  # 同次変換行列のリスト．0から順にT_Wo_BL, T_BL_0, T_0_1, ...
        
        self.T_Wo_j = []  # 同次変換行列のリスト．0から順にT_Wo_BL, T_Wo_0, T_Wo_1, ...
        for i, T in enumerate(self.T_i_j):
            if i == 0:
                self.T_Wo_j.append(T)
            else:
                self.T_Wo_j.append(self.T_Wo_j[i-1] @ T)
        
        return None
    
    
    def origins(self, q):
        """世界座標系から見た局所座標系の原点座標を返す"""
        
        # 同時変換行列を作成
        self.init_homo_transf_mat(q)
        origins = [np.zeros((3, 1))]
        for T in self.T_Wo_j:
            origins.append(T[0:3, 3:4])
        
        return origins
    
    

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D 
    
    kinema = MyCobotKinematics(
        d1 = 0.13156,
        a2 = 0.1104,
        a3 = 0.096,
        d4 = 0.06639,
        d5 = 0.07318,
        d6 = 0.0436
    )
    
    q = np.array([[pi/3, 0, 0, 0, 0, 0]]).T
    
    kinema.init_homo_transf_mat(q)
    x = kinema.origins(q)
    
    x_list = x[0]
    for temp in x[1:]:
        x_list = np.concatenate([x_list, temp], axis = 1)
    
    
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.plot(x_list[0, :], x_list[1, :], x_list[2, :])
    plt.show()