"""これが一番最新"""


import numpy as np
from math import cos, sin, tan, pi, sqrt



class BaxterKinematicsNew:
    """ローカル座標系の原点，ヤコビ行列等を計算\n
    ・行列べた書きを廃止
    """
    
    def __init__(self, **kwargs):
        self.L = kwargs.pop('L')
        self.h = kwargs.pop('h')
        self.H = kwargs.pop('H')
        self.L0 = kwargs.pop('L0')
        self.L1 = kwargs.pop('L1')
        self.L2 = kwargs.pop('L2')
        self.L3 = kwargs.pop('L3')
        self.L4 = kwargs.pop('L4')
        self.L5 = kwargs.pop('L5')
        self.L6 = kwargs.pop('L6')
    
    # def homo_transf_matrix(self, theta, a, alpha, d):
    #     """DH記法で使う
    #     ・未完成
    #     """
    #     return np.array([
    #         [cos(theta), -sin(theta), 0, a],
    #         [sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha)],
    #         [sin(theta)*sin(alpha), cos(theta)*sin(alpha), -cos(alpha), d*cos(alpha)],
    #         [0, 0, 0, 1],
    #     ])
    
    def init_homo_transf_mat(self, q):
        """同時変換行列を作成
        ・
        """
        q0, q1, q2, q3, q4, q5, q6, q7 = pi/4, q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        
        # 個別
        self.T_i_j = [
            np.array([
                [cos(q0), sin(q0), 0, self.L],
                [-sin(q0), cos(q0), 0, -self.h],
                [0, 0, 1, self.H],
                [0, 0, 0, 1],
            ]),
            np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, self.L0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q1), -sin(q1), 0, 0],
                [sin(q1), cos(q1), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [-sin(q2), -cos(q2), 0, self.L1],
                [0, 0, 1, 0],
                [-cos(q2), sin(q2), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q3), -sin(q3), 0, 0],
                [0, 0, -1, -self.L2],
                [sin(q3), cos(q3), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q4), -sin(q4), 0, self.L3],
                [0, 0, 1, 0],
                [-sin(q4), -cos(q4), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q5), -sin(q5), 0, 0],
                [0, 0, -1, -self.L4],
                [sin(q5), cos(q5), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q6), -sin(q6), 0, self.L5],
                [0, 0, 1, 0],
                [-sin(q6), -cos(q6), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [cos(q7), -sin(q7), 0, 0],
                [0, 0, -1, 0],
                [sin(q7), cos(q7), 0, 0],
                [0, 0, 0, 1],
            ]),
            np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, self.L6],
                [0, 0, 0, 1],
            ]),
        ]  # 同次変換行列のリスト．0から順にT_Wo_BL, T_BL_0, T_0_1, ...
        
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
    
    
    def init_jacobi(self, q):
        """ヤコビ行列を作成"""
        
        q0, q1, q2, q3, q4, q5, q6, q7 = pi/4, q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        
        dT_0_1dq1 = np.array([
            [-sin(q1), -cos(q1), 0, 0],
            [cos(q1), -sin(q1), 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ])
        dT_1_2dq2 = np.array([
            [-cos(q2), sin(q2), 0, 0],
            [0, 0, 0, 0],
            [sin(q2), cos(q2), 0, 0],
            [0, 0, 0, 0],
        ])
        dT_2_3dq3 = np.array([
            [-sin(q3), -cos(q3), 0, 0],
            [0, 0, 0, 0],
            [cos(q3), -sin(q3), 0, 0],
            [0, 0, 0, 0],
        ])
        dT_3_4dq4 = np.array([
            [-sin(q4), -cos(q4), 0, 0],
            [0, 0, 0, 0],
            [-cos(q4), sin(q4), 0, 0],
            [0, 0, 0, 0],
        ])
        dT_4_5dq5 = np.array([
            [-sin(q5), -cos(q5), 0, 0],
            [0, 0, 0, 0],
            [cos(q5), -sin(q5), 0, 0],
            [0, 0, 0, 0],
        ])
        dT_5_6dq6 = np.array([
            [-sin(q6), -cos(q6), 0, 0],
            [0, 0, 0, 0],
            [-cos(q6), sin(q6), 0, 0],
            [0, 0, 0, 0],
        ])
        dT_6_7dq7 = np.array([
            [-sin(q7), -cos(q7), 0, 0],
            [0, 0, 0, 0],
            [cos(q7), -sin(q7), 0, 0],
            [0, 0, 0, 0],
        ])
        
        dT_i_j = [
            dT_0_1dq1,
            dT_1_2dq2,
            dT_2_3dq3,
            dT_3_4dq4,
            dT_4_5dq5,
            dT_5_6dq6,
            dT_6_7dq7,
        ]
        #print(dT_i_j)
        
        dTdq1 = dTdq2 = dTdq3 = dTdq4 = dTdq5 = dTdq6 = dTdq7 = []
        for i, T in enumerate(self.T_i_j):
            if i == 0 or i == 1:
                dTdq1.append(np.zeros((4, 4)))
                dTdq2.append(np.zeros((4, 4)))
                dTdq3.append(np.zeros((4, 4)))
                dTdq4.append(np.zeros((4, 4)))
                dTdq5.append(np.zeros((4, 4)))
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
                
            elif i == 2:
                dTdq1.append(self.T_Wo_j[i-1] @ dT_0_1dq1)
                dTdq2.append(np.zeros((4, 4)))
                dTdq3.append(np.zeros((4, 4)))
                dTdq4.append(np.zeros((4, 4)))
                dTdq5.append(np.zeros((4, 4)))
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
            elif i == 3:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(self.T_Wo_j[i-1] @ dT_1_2dq2)
                dTdq3.append(np.zeros((4, 4)))
                dTdq4.append(np.zeros((4, 4)))
                dTdq5.append(np.zeros((4, 4)))
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
            elif i == 4:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(self.T_Wo_j[i-1] @ dT_2_3dq3)
                dTdq4.append(np.zeros((4, 4)))
                dTdq5.append(np.zeros((4, 4)))
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
            elif i == 5:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(dTdq3[i-1] @ T)
                dTdq4.append(self.T_Wo_j[i-1] @ dT_3_4dq4)
                dTdq5.append(np.zeros((4, 4)))
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
            elif i == 6:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(dTdq3[i-1] @ T)
                dTdq4.append(dTdq4[i-1] @ T)
                dTdq5.append(self.T_Wo_j[i-1] @ dT_4_5dq5)
                dTdq6.append(np.zeros((4, 4)))
                dTdq7.append(np.zeros((4, 4)))
            elif i == 7:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(dTdq3[i-1] @ T)
                dTdq4.append(dTdq4[i-1] @ T)
                dTdq5.append(dTdq5[i-1] @ T)
                dTdq6.append(self.T_Wo_j[i-1] @ dT_5_6dq6)
                dTdq7.append(np.zeros((4, 4)))
            elif i == 8:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(dTdq3[i-1] @ T)
                dTdq4.append(dTdq4[i-1] @ T)
                dTdq5.append(dTdq5[i-1] @ T)
                dTdq6.append(dTdq6[i-1] @ T)
                dTdq7.append(self.T_Wo_j[i-1] @ dT_6_7dq7)
            else:
                dTdq1.append(dTdq1[i-1] @ T)
                dTdq2.append(dTdq2[i-1] @ T)
                dTdq3.append(dTdq3[i-1] @ T)
                dTdq4.append(dTdq4[i-1] @ T)
                dTdq5.append(dTdq5[i-1] @ T)
                dTdq6.append(dTdq6[i-1] @ T)
                dTdq7.append(dTdq7[i-1] @ T)
            
            print('\n', i)
            print(dTdq1[-1], dTdq2[-1], dTdq3[-1], dTdq4[-1], dTdq5[-1], dTdq6[-1], dTdq7[-1])
        
        jacobi = [
            np.concatenate(
                [
                    dTdq1[i][0:3, 3:4],
                    dTdq2[i][0:3, 3:4],
                    dTdq3[i][0:3, 3:4],
                    dTdq4[i][0:3, 3:4],
                    dTdq5[i][0:3, 3:4],
                    dTdq6[i][0:3, 3:4],
                    dTdq7[i][0:3, 3:4],
                ],
                axis = 1,
            ) for i in range(11)
        ]
        self.jacobi = jacobi
        
        return None
    
    def calc_jacobi(self, q):
        self.init_jacobi(q)
        return self.jacobi


if __name__ == '__main__':
    BaxterKinema_2 = BaxterKinematicsNew(
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
    
    q1_init = 0
    q2_init = -31
    q3_init = 0
    q4_init = 43
    q5_init = 0
    q6_init = 72
    q7_init = 0
    q = np.array([[q1_init, q2_init, q3_init, q4_init, q5_init, q6_init, q7_init]]).T * np.pi / 180
    
    x_2 = BaxterKinema_2.origins(q)
    J_2 = BaxterKinema_2.calc_jacobi(q)
    #print(x_2)
    #print(J_2)