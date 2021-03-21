
import numpy as np
from math import pi, cos, sin, tan, atan


class OpemLabKineatics:
    """オープンラボのロボットアームのクラス"""
    def __init__(self, L0, L1, L2, L3, L4, L5, L6, H45):
        self.L0 = L0
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3 
        self.L4 = L4
        self.L5 = L5
        self.L6 = L6 
        self.H45 = H45
    
    # リンク基底原点位置
    def origin_0(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.zeros((3, 1))
        return z
    
    def origin_1(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [0], 
            [0], 
            [l0]
            ])
        return z
    
    def origin_2(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [0], 
            [0], 
            [l0 + l1]
            ])
        return z
    
    def origin_3(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [l2*sin(q2)*cos(q1)], 
            [l2*sin(q1)*sin(q2)], 
            [l0 + l1 + l2*cos(q2)]
            ])
        return z
    
    def origin_4(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [l2*sin(q2)*cos(q1) + l3*sin(q2 + q3)*cos(q1)], 
            [l2*sin(q1)*sin(q2) + l3*sin(q1)*sin(q2 + q3)], 
            [l0 + l1 + l2*cos(q2) + l3*cos(q2 + q3)]
            ])
        return z
    
    def origin_5(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [(l2*sin(q2) + l3*sin(q2 + q3) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*cos(q1)], 
            [(l2*sin(q2) + l3*sin(q2 + q3) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*sin(q1)], 
            [l0 + l1 + l2*cos(q2) + l3*cos(q2 + q3) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))]
            ])
        return z
    
    def origin_6(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [(l2*sin(q2) + l3*sin(q2 + q3) + l5*sin(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*cos(q1)], 
            [(l2*sin(q2) + l3*sin(q2 + q3) + l5*sin(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*sin(q1)], 
            [l0 + l1 + l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))]
            ])
        return z
    
    def origin_7(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [l2*sin(q2)*cos(q1) + l3*sin(q2 + q3)*cos(q1) + l5*sin(q2 + q3 + q4)*cos(q1) + l6*(sin(q1)*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3 + q4))*sin(q6) - l6*sin(q2 + q3 + q4)*cos(q1)*cos(q6) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1)], 
            [l2*sin(q1)*sin(q2) + l3*sin(q1)*sin(q2 + q3) + l5*sin(q1)*sin(q2 + q3 + q4) + l6*(sin(q1)*cos(q5)*cos(q2 + q3 + q4) - sin(q5)*cos(q1))*sin(q6) - l6*sin(q1)*sin(q2 + q3 + q4)*cos(q6) + (h45**2 + l4**2)**0.5*sin(q1)*sin(q2 + q3 + q4 + atan(h45/l4))], 
            [l0 + l1 + l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))]
            ])
        return z
    
    def origins(self, q):
        o_0 = self.origin_0(q)
        o_1 = self.origin_1(q)
        o_2 = self.origin_2(q)
        o_3 = self.origin_3(q)
        o_4 = self.origin_4(q)
        o_5 = self.origin_5(q)
        o_6 = self.origin_6(q)
        o_7 = self.origin_7(q)
        return [o_0, o_1, o_2, o_3, o_4, o_5, o_6, o_7]
    
    
    # 基底のヤコビ行列
    def jacobi_0(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def jacobi_1(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def jacobi_2(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def jacobi_3(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [
                -l2*sin(q1)*sin(q2), 
                l2*cos(q1)*cos(q2), 
                0, 
                0, 
                0, 
                0
            ], 
            [
                l2*sin(q2)*cos(q1), 
                l2*sin(q1)*cos(q2), 
                0, 
                0, 
                0, 
                0
            ], 
            [
                0, 
                -l2*sin(q2), 
                0, 
                0, 
                0, 
                0
            ]
            ])
        return z
    
    def jacobi_4(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [
                -(l2*sin(q2) + l3*sin(q2 + q3))*sin(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3))*cos(q1), 
                l3*cos(q1)*cos(q2 + q3), 
                0, 
                0, 
                0
            ], 
            [
                (l2*sin(q2) + l3*sin(q2 + q3))*cos(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3))*sin(q1), 
                l3*sin(q1)*cos(q2 + q3), 
                0, 
                0, 
                0
            ], 
            [
                0, 
                -l2*sin(q2) - l3*sin(q2 + q3), 
                -l3*sin(q2 + q3), 
                0, 
                0, 
                0
            ]
            ])
        return z
    
    def jacobi_5(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [
                -(l2*sin(q2) + l3*sin(q2 + q3) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l3*cos(q2 + q3) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (h45**2 + l4**2)**0.5*cos(q1)*cos(q2 + q3 + q4 + atan(h45/l4)), 
                0, 
                0
            ], 
            [
                (l2*sin(q2) + l3*sin(q2 + q3) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l3*cos(q2 + q3) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)), 
                0, 
                0
            ], 
            [
                0, 
                -l2*sin(q2) - l3*sin(q2 + q3) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -l3*sin(q2 + q3) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                0, 
                0
            ]
            ])
        return z
    
    def jacobi_6(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [
                -(l2*sin(q2) + l3*sin(q2 + q3) + l5*sin(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                0, 
                0
            ], 
            [
                (l2*sin(q2) + l3*sin(q2 + q3) + l5*sin(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                0, 
                0
            ], 
            [
                0, 
                -l2*sin(q2) - l3*sin(q2 + q3) - l5*sin(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -l3*sin(q2 + q3) - l5*sin(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -l5*sin(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                0, 
                0
            ]
            ])
        return z
    
    def jacobi_7(self, q):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        z = np.array([
            [
                -l2*sin(q1)*sin(q2) - l3*sin(q1)*sin(q2 + q3) - l5*sin(q1)*sin(q2 + q3 + q4) - l6*(sin(q1)*cos(q5)*cos(q2 + q3 + q4) - sin(q5)*cos(q1))*sin(q6) + l6*sin(q1)*sin(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*sin(q1)*sin(q2 + q3 + q4 + atan(h45/l4)), 
                (l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                (l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*cos(q1), 
                l6*(sin(q1)*cos(q5) - sin(q5)*cos(q1)*cos(q2 + q3 + q4))*sin(q6), 
                l6*((sin(q1)*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3 + q4))*cos(q6) + sin(q6)*sin(q2 + q3 + q4)*cos(q1))
            ], 
            [
                l2*sin(q2)*cos(q1) + l3*sin(q2 + q3)*cos(q1) + l5*sin(q2 + q3 + q4)*cos(q1) + l6*(sin(q1)*sin(q5) + cos(q1)*cos(q5)*cos(q2 + q3 + q4))*sin(q6) - l6*sin(q2 + q3 + q4)*cos(q1)*cos(q6) + (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), 
                (l2*cos(q2) + l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l3*cos(q2 + q3) + l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                (l5*cos(q2 + q3 + q4) - l6*sin(q6)*sin(q2 + q3 + q4)*cos(q5) - l6*cos(q6)*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*sin(q1), 
                -l6*(sin(q1)*sin(q5)*cos(q2 + q3 + q4) + cos(q1)*cos(q5))*sin(q6), 
                l6*((sin(q1)*cos(q5)*cos(q2 + q3 + q4) - sin(q5)*cos(q1))*cos(q6) + sin(q1)*sin(q6)*sin(q2 + q3 + q4))
            ], 
            [
                0, 
                -l2*sin(q2) - l3*sin(q2 + q3) - l5*sin(q2 + q3 + q4) - l6*sin(q6)*cos(q5)*cos(q2 + q3 + q4) + l6*sin(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -l3*sin(q2 + q3) - l5*sin(q2 + q3 + q4) - l6*sin(q6)*cos(q5)*cos(q2 + q3 + q4) + l6*sin(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), 
                -l5*sin(q2 + q3 + q4) - l6*sin(q6)*cos(q5)*cos(q2 + q3 + q4) + l6*sin(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4)), l6*sin(q5)*sin(q6)*sin(q2 + q3 + q4), 
                l6*(sin(q6)*cos(q2 + q3 + q4) - sin(q2 + q3 + q4)*cos(q5)*cos(q6))
            ]
            ])
        return z
    
    def jacobi_all(self, q):
        j0 = self.jacobi_0(q)
        j1 = self.jacobi_1(q)
        j2 = self.jacobi_2(q)
        j3 = self.jacobi_3(q)
        j4 = self.jacobi_4(q)
        j5 = self.jacobi_5(q)
        j6 = self.jacobi_6(q)
        j7 = self.jacobi_7(q)
        return [j0, j1, j2, j3, j4, j5, j6, j7]
    
    
    # 基底のヤコビ行列の時間微分
    def djacobi_0(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def djacobi_1(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def djacobi_2(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.zeros((3, 6))
        return z
    
    def djacobi_3(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.array([
            [-l2*(dq1*sin(q2)*cos(q1) + dq2*sin(q1)*cos(q2)), -l2*(dq1*sin(q1)*cos(q2) + dq2*sin(q2)*cos(q1)), 0, 0, 0, 0], 
            [l2*(-dq1*sin(q1)*sin(q2) + dq2*cos(q1)*cos(q2)), l2*(dq1*cos(q1)*cos(q2) - dq2*sin(q1)*sin(q2)), 0, 0, 0, 0], 
            [0, -dq2*l2*cos(q2), 0, 0, 0, 0]
            ])
        return z
    
    def djacobi_4(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.array([
            [-dq1*l2*sin(q2)*cos(q1) - dq1*l3*sin(q2 + q3)*cos(q1) - dq2*l2*sin(q1)*cos(q2) - l3*(dq2 + dq3)*sin(q1)*cos(q2 + q3), -dq1*l2*sin(q1)*cos(q2) - dq1*l3*sin(q1)*cos(q2 + q3) - dq2*l2*sin(q2)*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1), -l3*(dq1*sin(q1)*cos(q2 + q3) + (dq2 + dq3)*sin(q2 + q3)*cos(q1)), 0, 0, 0], 
            [-dq1*l2*sin(q1)*sin(q2) - dq1*l3*sin(q2 + q3)*sin(q1) + dq2*l2*cos(q1)*cos(q2) + l3*(dq2 + dq3)*cos(q2 + q3)*cos(q1), dq1*l2*cos(q1)*cos(q2) + dq1*l3*cos(q2 + q3)*cos(q1) - dq2*l2*sin(q1)*sin(q2) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1), l3*(dq1*cos(q2 + q3)*cos(q1) - (dq2 + dq3)*sin(q2 + q3)*sin(q1)), 0, 0, 0], 
            [0, -dq2*l2*cos(q2) - l3*(dq2 + dq3)*cos(q2 + q3), -l3*(dq2 + dq3)*cos(q2 + q3), 0, 0, 0]
            ])
        return z
    
    def djacobi_5(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.array([
            [-dq1*l2*sin(q2)*cos(q1) - dq1*l3*sin(q2 + q3)*cos(q1) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*cos(q2) - l3*(dq2 + dq3)*sin(q1)*cos(q2 + q3) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)), -dq1*l2*sin(q1)*cos(q2) - dq1*l3*sin(q1)*cos(q2 + q3) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - dq2*l2*sin(q2)*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -dq1*l3*sin(q1)*cos(q2 + q3) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -(h45**2 + l4**2)**0.5*(dq1*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1)), 0, 0], 
            [-dq1*l2*sin(q1)*sin(q2) - dq1*l3*sin(q2 + q3)*sin(q1) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1) + dq2*l2*cos(q1)*cos(q2) + l3*(dq2 + dq3)*cos(q2 + q3)*cos(q1) + (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1), dq1*l2*cos(q1)*cos(q2) + dq1*l3*cos(q2 + q3)*cos(q1) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*sin(q2) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), dq1*l3*cos(q2 + q3)*cos(q1) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), (h45**2 + l4**2)**0.5*(dq1*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - (dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1)), 0, 0], 
            [0, -dq2*l2*cos(q2) - l3*(dq2 + dq3)*cos(q2 + q3) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), -l3*(dq2 + dq3)*cos(q2 + q3) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), -(h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), 0, 0]
            ])
        return z
    
    def djacobi_6(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.array([
            [-dq1*l2*sin(q2)*cos(q1) - dq1*l3*sin(q2 + q3)*cos(q1) - dq1*l5*sin(q2 + q3 + q4)*cos(q1) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*cos(q2) - l3*(dq2 + dq3)*sin(q1)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)), -dq1*l2*sin(q1)*cos(q2) - dq1*l3*sin(q1)*cos(q2 + q3) - dq1*l5*sin(q1)*cos(q2 + q3 + q4) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - dq2*l2*sin(q2)*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -dq1*l3*sin(q1)*cos(q2 + q3) - dq1*l5*sin(q1)*cos(q2 + q3 + q4) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -dq1*l5*sin(q1)*cos(q2 + q3 + q4) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), 0, 0], 
            [-dq1*l2*sin(q1)*sin(q2) - dq1*l3*sin(q2 + q3)*sin(q1) - dq1*l5*sin(q2 + q3 + q4)*sin(q1) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1) + dq2*l2*cos(q1)*cos(q2) + l3*(dq2 + dq3)*cos(q2 + q3)*cos(q1) + l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q1) + (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1), dq1*l2*cos(q1)*cos(q2) + dq1*l3*cos(q2 + q3)*cos(q1) + dq1*l5*cos(q2 + q3 + q4)*cos(q1) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*sin(q2) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), dq1*l3*cos(q2 + q3)*cos(q1) + dq1*l5*cos(q2 + q3 + q4)*cos(q1) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), dq1*l5*cos(q2 + q3 + q4)*cos(q1) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), 0, 0], 
            [0, -dq2*l2*cos(q2) - l3*(dq2 + dq3)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), -l3*(dq2 + dq3)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), -(l5*cos(q2 + q3 + q4) + (h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4)))*(dq2 + dq3 + dq4), 0, 0]])
        return z
    
    def djacobi_7(self, q, dq):
        l0, l1, l2, l3, l4, l5, l6, h45 = self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6, self.H45
        q1, q2, q3, q4, q5, q6 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0]
        dq1, dq2, dq3, dq4, dq5, dq6 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0]
        z = np.array([
            [-dq1*l2*sin(q2)*cos(q1) - dq1*l3*sin(q2 + q3)*cos(q1) - dq1*l5*sin(q2 + q3 + q4)*cos(q1) + dq1*l6*sin(q2 + q3 + q4)*cos(q1)*cos(q6) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*cos(q2) - dq6*l6*(sin(q1)*cos(q2 + q3 + q4)*cos(q5) - sin(q5)*cos(q1))*cos(q6) - dq6*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q6) - l3*(dq2 + dq3)*sin(q1)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4) + l6*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4)*cos(q6) + l6*(-dq1*sin(q1)*sin(q5) - dq1*cos(q2 + q3 + q4)*cos(q1)*cos(q5) + dq5*sin(q1)*sin(q5)*cos(q2 + q3 + q4) + dq5*cos(q1)*cos(q5) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*cos(q5))*sin(q6) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)), -dq1*l2*sin(q1)*cos(q2) - dq1*l3*sin(q1)*cos(q2 + q3) - dq1*l5*sin(q1)*cos(q2 + q3 + q4) + dq1*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q6)*cos(q5) + dq1*l6*sin(q1)*cos(q2 + q3 + q4)*cos(q6) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) - dq2*l2*sin(q2)*cos(q1) + dq5*l6*sin(q2 + q3 + q4)*sin(q5)*sin(q6)*cos(q1) - dq6*l6*sin(q2 + q3 + q4)*cos(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q6)*cos(q2 + q3 + q4)*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q6)*cos(q2 + q3 + q4)*cos(q1)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -dq1*l3*sin(q1)*cos(q2 + q3) - dq1*l5*sin(q1)*cos(q2 + q3 + q4) + dq1*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q6)*cos(q5) + dq1*l6*sin(q1)*cos(q2 + q3 + q4)*cos(q6) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) + dq5*l6*sin(q2 + q3 + q4)*sin(q5)*sin(q6)*cos(q1) - dq6*l6*sin(q2 + q3 + q4)*cos(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q6)*cos(q2 + q3 + q4)*cos(q1) - l3*(dq2 + dq3)*sin(q2 + q3)*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q6)*cos(q2 + q3 + q4)*cos(q1)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), -dq1*l5*sin(q1)*cos(q2 + q3 + q4) + dq1*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q6)*cos(q5) + dq1*l6*sin(q1)*cos(q2 + q3 + q4)*cos(q6) - dq1*(h45**2 + l4**2)**0.5*sin(q1)*cos(q2 + q3 + q4 + atan(h45/l4)) + dq5*l6*sin(q2 + q3 + q4)*sin(q5)*sin(q6)*cos(q1) - dq6*l6*sin(q2 + q3 + q4)*cos(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q6)*cos(q2 + q3 + q4)*cos(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q6)*cos(q2 + q3 + q4)*cos(q1)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*cos(q1), l6*(dq6*(sin(q1)*cos(q5) - sin(q5)*cos(q2 + q3 + q4)*cos(q1))*cos(q6) + (dq1*sin(q1)*sin(q5)*cos(q2 + q3 + q4) + dq1*cos(q1)*cos(q5) - dq5*sin(q1)*sin(q5) - dq5*cos(q2 + q3 + q4)*cos(q1)*cos(q5) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q5)*cos(q1))*sin(q6)), l6*(-dq1*sin(q2 + q3 + q4)*sin(q1)*sin(q6) - dq6*(sin(q1)*sin(q5) + cos(q2 + q3 + q4)*cos(q1)*cos(q5))*sin(q6) + dq6*sin(q2 + q3 + q4)*cos(q1)*cos(q6) + (dq2 + dq3 + dq4)*sin(q6)*cos(q2 + q3 + q4)*cos(q1) - (dq1*sin(q1)*cos(q2 + q3 + q4)*cos(q5) - dq1*sin(q5)*cos(q1) - dq5*sin(q1)*cos(q5) + dq5*sin(q5)*cos(q2 + q3 + q4)*cos(q1) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1)*cos(q5))*cos(q6))], 
            [-dq1*l2*sin(q1)*sin(q2) - dq1*l3*sin(q2 + q3)*sin(q1) - dq1*l5*sin(q2 + q3 + q4)*sin(q1) + dq1*l6*sin(q2 + q3 + q4)*sin(q1)*cos(q6) - dq1*(h45**2 + l4**2)**0.5*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1) + dq2*l2*cos(q1)*cos(q2) + dq6*l6*(sin(q1)*sin(q5) + cos(q2 + q3 + q4)*cos(q1)*cos(q5))*cos(q6) + dq6*l6*sin(q2 + q3 + q4)*sin(q6)*cos(q1) + l3*(dq2 + dq3)*cos(q2 + q3)*cos(q1) + l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q1) - l6*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q1)*cos(q6) - l6*(dq1*sin(q1)*cos(q2 + q3 + q4)*cos(q5) - dq1*sin(q5)*cos(q1) - dq5*sin(q1)*cos(q5) + dq5*sin(q5)*cos(q2 + q3 + q4)*cos(q1) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*cos(q1)*cos(q5))*sin(q6) + (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1), dq1*l2*cos(q1)*cos(q2) + dq1*l3*cos(q2 + q3)*cos(q1) + dq1*l5*cos(q2 + q3 + q4)*cos(q1) - dq1*l6*sin(q2 + q3 + q4)*sin(q6)*cos(q1)*cos(q5) - dq1*l6*cos(q2 + q3 + q4)*cos(q1)*cos(q6) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) - dq2*l2*sin(q1)*sin(q2) + dq5*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q5)*sin(q6) - dq6*l6*sin(q2 + q3 + q4)*sin(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q1)*sin(q6)*cos(q2 + q3 + q4) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q1)*sin(q6)*cos(q2 + q3 + q4)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), dq1*l3*cos(q2 + q3)*cos(q1) + dq1*l5*cos(q2 + q3 + q4)*cos(q1) - dq1*l6*sin(q2 + q3 + q4)*sin(q6)*cos(q1)*cos(q5) - dq1*l6*cos(q2 + q3 + q4)*cos(q1)*cos(q6) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) + dq5*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q5)*sin(q6) - dq6*l6*sin(q2 + q3 + q4)*sin(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q1)*sin(q6)*cos(q2 + q3 + q4) - l3*(dq2 + dq3)*sin(q2 + q3)*sin(q1) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q1)*sin(q6)*cos(q2 + q3 + q4)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), dq1*l5*cos(q2 + q3 + q4)*cos(q1) - dq1*l6*sin(q2 + q3 + q4)*sin(q6)*cos(q1)*cos(q5) - dq1*l6*cos(q2 + q3 + q4)*cos(q1)*cos(q6) + dq1*(h45**2 + l4**2)**0.5*cos(q2 + q3 + q4 + atan(h45/l4))*cos(q1) + dq5*l6*sin(q2 + q3 + q4)*sin(q1)*sin(q5)*sin(q6) - dq6*l6*sin(q2 + q3 + q4)*sin(q1)*cos(q5)*cos(q6) + dq6*l6*sin(q1)*sin(q6)*cos(q2 + q3 + q4) - l5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*cos(q6) - l6*(dq2 + dq3 + dq4)*sin(q1)*sin(q6)*cos(q2 + q3 + q4)*cos(q5) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4 + atan(h45/l4))*sin(q1), l6*(-dq6*(sin(q1)*sin(q5)*cos(q2 + q3 + q4) + cos(q1)*cos(q5))*cos(q6) + (dq1*sin(q1)*cos(q5) - dq1*sin(q5)*cos(q2 + q3 + q4)*cos(q1) - dq5*sin(q1)*cos(q2 + q3 + q4)*cos(q5) + dq5*sin(q5)*cos(q1) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*sin(q5))*sin(q6)), l6*(dq1*sin(q2 + q3 + q4)*sin(q6)*cos(q1) - dq6*(sin(q1)*cos(q2 + q3 + q4)*cos(q5) - sin(q5)*cos(q1))*sin(q6) + dq6*sin(q2 + q3 + q4)*sin(q1)*cos(q6) + (dq2 + dq3 + dq4)*sin(q1)*sin(q6)*cos(q2 + q3 + q4) - (-dq1*sin(q1)*sin(q5) - dq1*cos(q2 + q3 + q4)*cos(q1)*cos(q5) + dq5*sin(q1)*sin(q5)*cos(q2 + q3 + q4) + dq5*cos(q1)*cos(q5) + (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q1)*cos(q5))*cos(q6))], 
            [0, -dq2*l2*cos(q2) + dq5*l6*sin(q5)*sin(q6)*cos(q2 + q3 + q4) - dq6*l6*sin(q2 + q3 + q4)*sin(q6) - dq6*l6*cos(q2 + q3 + q4)*cos(q5)*cos(q6) - l3*(dq2 + dq3)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q6)*cos(q5) + l6*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), dq5*l6*sin(q5)*sin(q6)*cos(q2 + q3 + q4) - dq6*l6*sin(q2 + q3 + q4)*sin(q6) - dq6*l6*cos(q2 + q3 + q4)*cos(q5)*cos(q6) - l3*(dq2 + dq3)*cos(q2 + q3) - l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q6)*cos(q5) + l6*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), dq5*l6*sin(q5)*sin(q6)*cos(q2 + q3 + q4) - dq6*l6*sin(q2 + q3 + q4)*sin(q6) - dq6*l6*cos(q2 + q3 + q4)*cos(q5)*cos(q6) - l5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4) + l6*(dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q6)*cos(q5) + l6*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q6) - (h45**2 + l4**2)**0.5*(dq2 + dq3 + dq4)*cos(q2 + q3 + q4 + atan(h45/l4)), l6*(dq5*sin(q2 + q3 + q4)*sin(q6)*cos(q5) + dq6*sin(q2 + q3 + q4)*sin(q5)*cos(q6) + (dq2 + dq3 + dq4)*sin(q5)*sin(q6)*cos(q2 + q3 + q4)), l6*(dq5*sin(q2 + q3 + q4)*sin(q5)*cos(q6) + dq6*sin(q2 + q3 + q4)*sin(q6)*cos(q5) + dq6*cos(q2 + q3 + q4)*cos(q6) - (dq2 + dq3 + dq4)*sin(q2 + q3 + q4)*sin(q6) - (dq2 + dq3 + dq4)*cos(q2 + q3 + q4)*cos(q5)*cos(q6))]])
        return z
    
    def djacobi_all(self, q, dq):
        """ヤコビ行列の時間微分を全部返す"""
        dj0 = self.djacobi_0(q, dq)
        dj1 = self.djacobi_1(q, dq)
        dj2 = self.djacobi_2(q, dq)
        dj3 = self.djacobi_3(q, dq)
        dj4 = self.djacobi_4(q, dq)
        dj5 = self.djacobi_5(q, dq)
        dj6 = self.djacobi_6(q, dq)
        dj7 = self.djacobi_7(q, dq)
        return [dj0, dj1, dj2, dj3, dj4, dj5, dj6, dj7]


