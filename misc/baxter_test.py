import numpy as np
from math import cos, sin, tan, pi, sqrt
import time


class Kinematics:
    """ローカル座標系の原点，ヤコビ行列等を計算
    ・ちょっとはやい
    """
    
    L = 278e-3
    h = 64e-3
    H = 1104e-3
    L0 = 270.35e-3
    L1 = 69e-3
    L2 = 364.35e-3
    L3 = 69e-3
    L4 = 374.29e-3
    L5 = 10e-3
    L6 = 368.3e-3
    
    q_min_ = [-143, -123, -173,   -3, -175,  -90, -175]
    q_max_ = [  51,   60,  173,  150,  175,  120,  175]
    q_min = [q * pi / 180 for q in q_min_]
    q_max = [q * pi / 180 for q in q_max_]
    
    def __init__(self, **kwargs):
        """"""
        return
    
    
    def homo_transf_matrix(self, alpha, a, d, theta):
        """同次変換行列(DH記法)
        """
        return np.array([
            [cos(theta), -sin(theta), 0, a],
            [sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha)],
            [sin(theta)*sin(alpha), cos(theta)*sin(alpha), -cos(alpha), d*cos(alpha)],
            [0, 0, 0, 1],
            ])
    
    
    def diff_homo_transf_matrix_by_theta(self, alpha, a, d, theta):
        """同次変換行列をthetaで微分"""
        return np.array([
            [-sin(theta), -cos(theta), 0, 0],
            [cos(theta)*cos(alpha), -sin(theta)*cos(alpha), 0, 0],
            [cos(theta)*sin(alpha), -sin(theta)*sin(alpha), 0, 0],
            [0, 0, 0, 0],
            ])
    
    
    def update_homo_transf_matrix(self, q):
        """同時変換行列を更新
        ・
        """
        q = np.ravel(q).tolist()
        
        DH_params = [
            (    0,       0, self.L0,         0),  # T_BL_0
            (    0,       0,       0,      q[0]),
            (-pi/2, self.L1,       0, q[1]+pi/2),
            ( pi/2,       0, self.L2,      q[2]),
            (-pi/2, self.L3,       0,      q[3]),
            ( pi/2,       0, self.L4,      q[4]),
            (-pi/2, self.L5,       0,      q[5]),
            ( pi/2,       0,       0,      q[6]),
        ]
        
        def T_Wo_BL():
            """これだけDHパラメータわからん"""
            return np.array([
                [cos(pi/4), sin(pi/4), 0, self.L],
                [-sin(pi/4), cos(pi/4), 0, -self.h],
                [0, 0, 1, self.H],
                [0, 0, 0, 1],
            ])
        
        # 同次変換行列のリスト．0から順にT_Wo_BL, T_BL_0, T_0_1, ...
        self.T_i_j = [T_Wo_BL()] + [self.homo_transf_matrix(*param) for param in DH_params]
        
        self.T_Wo_j = []  # 同次変換行列のリスト．0から順にT_Wo_BL, T_Wo_0, T_Wo_1, ...
        for i, T in enumerate(self.T_i_j):
            if i == 0:
                self.T_Wo_j.append(T)
            else:
                self.T_Wo_j.append(self.T_Wo_j[i-1] @ T)
        
        return
    
    
    def update_diff_homo_transf_mat_by_theta(self, q):
        """同次変換行列のtheta微分を更新"""
        
        q = np.ravel(q).tolist()
        
        DH_params = [
            (    0,       0, self.L0,         0),  # T_BL_0
            (    0,       0,       0,      q[0]),
            (-pi/2, self.L1,       0, q[1]+pi/2),
            ( pi/2,       0, self.L2,      q[2]),
            (-pi/2, self.L3,       0,      q[3]),
            ( pi/2,       0, self.L4,      q[4]),
            (-pi/2, self.L5,       0,      q[5]),
            ( pi/2,       0,       0,      q[6]),
        ]
        
        def dT_Wo_BL_dq():
            """これだけDHパラメータわからん"""
            return np.zeros((4, 4))
        
        # 同次変換行列のtheta微分のリスト．0から順にdT_Wo_BL, dT_BL_0, dT_0_1, ...
        self.dT_i_j_dq = [dT_Wo_BL_dq()] + [self.diff_homo_transf_matrix_by_theta(*param) for param in DH_params]
        
        dTdq1 = [self.dT_i_j_dq[0] @ self.T_Wo_j[1]]
        for T in self.T_i_j[3:]:
            dTdq1.append(T @ dTdq1[-1])
        
        dTdq2 = [np.zeros((4, 4))]
        dTdq2.append(self.dT_i_j_dq[1] @ self.T_i_j[2])
        for T in self.T_i_j[4:]:
            dTdq2.append(T @ dTdq2[-1])
        
        dTdq3 = [np.zeros((4, 4))] * 2
        dTdq3.append(self.dT_i_j_dq[2] @ self.T_i_j[3])
        for T in self.T_i_j[5:]:
            dTdq3.append(T @ dTdq3[-1])
        
        dTdq4 = [np.zeros((4, 4))] * 3
        dTdq4.append(self.dT_i_j_dq[3] @ self.T_i_j[4])
        for T in self.T_i_j[6:]:
            dTdq4.append(T @ dTdq4[-1])
        
        dTdq5 = [np.zeros((4, 4))] * 4
        dTdq5.append(self.dT_i_j_dq[4] @ self.T_i_j[5])
        for T in self.T_i_j[7:]:
            dTdq5.append(T @ dTdq5[-1])
        
        dTdq6 = [np.zeros((4, 4))] * 5
        dTdq6.append(self.dT_i_j_dq[5] @ self.T_i_j[6])
        for T in self.T_i_j[8:]:
            dTdq6.append(T @ dTdq6[-1])
        
        dTdq7 = [np.zeros((4, 4))] * 6
        dTdq7.append(self.dT_i_j_dq[6] @ self.T_i_j[7])
        for T in self.T_i_j[9:]:
            dTdq7.append(T @ dTdq7[-1])
        
        dTdq = [dTdq1, dTdq2, dTdq3, dTdq4, dTdq5, dTdq6, dTdq7,]
        
        print(len(dTdq[1]))
        
        jacobi = []
        for i in range(len(dTdq1)):
            J = []
            for j in range(7):
                J.append(dTdq[j][i][0:3, 3:4])
            jacobi.append(np.concatenate(J, axis = 1))
        self.jacobi = jacobi
        
        return
        
    
    
    
    
    
    def o(self,):
        """基底"""
        
        o_0 = [np.zeros((3, 1))]
        o_1_ = [T[0:3, 3:4] for T in self.T_Wo_j]
        os = o_0 + o_1_
        
        return os