"""これが一番最新"""


import numpy as np
from math import cos, sin, tan, pi, sqrt



class BaxterKinematicsNew:
    """ローカル座標系の原点，ヤコビ行列等を計算
    ・
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
        q0, q1, q2, q3, q4, q5, q6, q7 = pi/4, q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        
        # 個別
        self.T_Wo_BL = np.array([
            [cos(q0), sin(q0), 0, self.L],
            [-sin(q0), cos(q0), 0, -self.h],
            [0, 0, 1, self.H],
            [0, 0, 0, 1],
        ])
        self.T_BL_0 = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, self.L0],
            [0, 0, 0, 1],
        ])
        self.T_0_1 = np.array([
            [cos(q1), -sin(q1), 0, 0],
            [sin(q1), cos(q1), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ])
        self.T_1_2 = np.array([
            [-sin(q2), -cos(q2), 0, self.L1],
            [0, 0, 1, 0],
            [-cos(q2), sin(q2), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_2_3 = np.array([
            [cos(q3), -sin(q3), 0, 0],
            [0, 0, -1, -self.L2],
            [sin(q3), cos(q3), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_3_4 = np.array([
            [cos(q4), -sin(q4), 0, self.L3],
            [0, 0, 1, 0],
            [-sin(q4), -cos(q4), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_4_5 = np.array([
            [cos(q5), -sin(q5), 0, 0],
            [0, 0, -1, -self.L4],
            [sin(q5), cos(q5), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_5_6 = np.array(
            [[cos(q6), -sin(q6), 0, self.L5],
            [0, 0, 1, 0],
            [-sin(q6), -cos(q6), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_6_7 = np.array([
            [cos(q7), -sin(q7), 0, 0],
            [0, 0, -1, 0],
            [sin(q7), cos(q7), 0, 0],
            [0, 0, 0, 1],
        ])
        self.T_7_GL = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, self.L6],
            [0, 0, 0, 1],
        ])
        
        # 全体
        self.T_Wo_0 = self.T_Wo_BL @ self.T_BL_0
        self.T_Wo_1 = self.T_Wo_0 @ self.T_0_1
        self.T_Wo_2 = self.T_Wo_1 @ self.T_1_2
        self.T_Wo_3 = self.T_Wo_2 @ self.T_2_3
        self.T_Wo_4 = self.T_Wo_3 @ self.T_3_4
        self.T_Wo_5 = self.T_Wo_4 @ self.T_4_5
        self.T_Wo_6 = self.T_Wo_5 @ self.T_5_6
        self.T_Wo_7 = self.T_Wo_6 @ self.T_6_7
        self.T_Wo_GL = self.T_Wo_7 @ self.T_7_GL
    
    
    def calc_origins(self, q):
        """世界座標系から見た局所座標系の原点座標を返す"""
        
        # 同時変換行列を作成
        self.init_homo_transf_mat(q)
        
        return [
            np.zeros((3, 1)),
            self.T_Wo_BL[0:3, 3:4],
            self.T_Wo_0[0:3, 3:4],
            self.T_Wo_1[0:3, 3:4],
            self.T_Wo_2[0:3, 3:4],
            self.T_Wo_3[0:3, 3:4],
            self.T_Wo_4[0:3, 3:4],
            self.T_Wo_5[0:3, 3:4],
            self.T_Wo_6[0:3, 3:4],
            self.T_Wo_7[0:3, 3:4],
            self.T_Wo_GL[0:3, 3:4],
        ]

        # 計算

        
        
        
        
        
        # # 計算
        # homo_transf_mat_all = [
        #     self.T_BL_0,
        #     self.T_0_1,
        #     self.T_1_2,
        #     self.T_2_3,
        #     self.T_3_4,
        #     self.T_4_5,
        #     self.T_5_6,
        #     self.T_6_7,
        #     self.T_7_GL,
        # ]
        
        # name_temp = []
        # temp = self.T_Wo_BL
        # for mat in homo_transf_mat_all:
        #     temp = temp @ mat
        #     self.name_temp[]