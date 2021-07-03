"""baxterのタスク写像"""


import numpy as np

import baxter_oldold


class Maps(baxter_oldold.BaxterFuncs):
    """baxterの写像"""
    
    def __init__(self,):
        
        super().__init__()
        self.init_r_bars()
        
        self.func_T = [
            self.T0, self.T1, self.T2, self.T3, self.T4,
            self.T5, self.T6, self.T7, self.T8, self.T9,
        ]  # 同次変換行列（関数）のリスト
        
        self.func_Jax = [
            self.Jax0, self.Jax1, self.Jax2, self.Jax3, self.Jax4,
            self.Jax5, self.Jax6, self.Jax7, self.Jax8, self.Jax9,
        ]
        
        self.func_Jay = [
            self.Jay0, self.Jay1, self.Jay2, self.Jay3, self.Jay4,
            self.Jay5, self.Jay6, self.Jay7, self.Jay8, self.Jay9,
        ]
        
        self.func_Jaz = [
            self.Jaz0, self.Jaz1, self.Jaz2, self.Jaz3, self.Jaz4,
            self.Jaz5, self.Jaz6, self.Jaz7, self.Jaz8, self.Jaz9,
        ]
        
        self.func_Jo = [
            self.Jo0, self.Jo1, self.Jo2, self.Jo3, self.Jo4,
            self.Jo5, self.Jo6, self.Jo7, self.Jo8, self.Jo9,
        ]
        
        self.func_dJax = [
            self.dJax0, self.dJax1, self.dJax2, self.dJax3, self.dJax4,
            self.dJax5, self.dJax6, self.dJax7, self.dJax8, self.dJax9,
        ]
        
        self.func_dJay = [
            self.dJay0, self.dJay1, self.dJay2, self.dJay3, self.dJay4,
            self.dJay5, self.dJay6, self.dJay7, self.dJay8, self.dJay9,
        ]
        
        self.func_dJaz = [
            self.dJaz0, self.dJaz1, self.dJaz2, self.dJaz3, self.dJaz4,
            self.dJaz5, self.dJaz6, self.dJaz7, self.dJaz8, self.dJaz9,
        ]
        
        self.func_dJo = [
            self.dJo0, self.dJo1, self.dJo2, self.dJo3, self.dJo4,
            self.dJo5, self.dJo6, self.dJo7, self.dJo8, self.dJo9,
        ]
        
        
        self.maps = {}
    
    
    def init_r_bars(self,):
        """制御点座標を追加"""
        
        r_bar_1 = [
            np.array([[0, self.L1/2, -self.L0/2]]).T,
            np.array([[0, -self.L1/2, -self.L0/2]]).T,
            np.array([[self.L1/2, 0, -self.L0/2]]).T,
            np.array([[-self.L1/2, 0, -self.L0/2]]).T,
        ]  # 1座標系からみた制御点位置
        
        r_bar_2 = [
            np.array([[0, 0, self.L3/2]]).T,
            np.array([[0, 0, -self.L3/2]]).T,
        ]
        
        r_bar_3 = [
            np.array([[0, self.L3/2, -self.L2*2/3]]).T,
            np.array([[0, -self.L3/2, -self.L2*2/3]]).T,
            np.array([[self.L3/2, 0, -self.L2*2/3]]).T,
            np.array([[-self.L3/2, 0, -self.L2*2/3]]).T,
            np.array([[0, self.L3/2, -self.L2*1/3]]).T,
            np.array([[0, -self.L3/2, -self.L2*1/3]]).T,
            np.array([[self.L3/2, 0, -self.L2*1/3]]).T,
            np.array([[-self.L3/2, 0, -self.L2*1/3]]).T,
        ]
        
        r_bar_4 = [
            np.array([[0, 0, self.L3/2]]).T,
            np.array([[0, 0, -self.L3/2]]).T,
        ]
        
        r_bar_5 = [
            np.array([[0, self.L5/2, -self.L4/2]]).T,
            np.array([[0, -self.L5/2, -self.L4/2]]).T,
            np.array([[self.L5/2, 0, -self.L4/2]]).T,
            np.array([[-self.L5/2, 0, -self.L4/2]]).T,
        ]
        
        r_bar_6 = [
            np.array([[0, 0, self.L5/2]]).T,
            np.array([[0, 0, -self.L5/2]]).T,
        ]
        
        r_bar_7 = [
            np.array([[0, self.L5/2, self.L6/2]]).T,
            np.array([[0, -self.L5/2, self.L6/2]]).T,
            np.array([[self.L5/2, 0, self.L6/2]]).T,
            np.array([[-self.L5/2, 0, self.L6/2]]).T,
        ]
        
        r_bar_GL = [
            np.array([[0, 0, 0]]).T
        ]
        
        # 追加
        self.r_bars = [
            r_bar_1, r_bar_2, r_bar_3, r_bar_4, r_bar_5, r_bar_6, r_bar_7, r_bar_GL,
        ]
        return
    
    def update_all_baxter(self,):
        self.T = [f(self.q_list) for f in self.func_T]
        self.Jax = [f(self.q_list) for f in self.func_Jax]
        self.Jay = [f(self.q_list) for f in self.func_Jay]
        self.Jaz = [f(self.q_list) for f in self.func_Jaz]
        self.Jo = [f(self.q_list) for f in self.func_Jo]
        self.dJax = [f(self.q_list, self.dq_list) for f in self.func_dJax]
        self.dJay = [f(self.q_list, self.dq_list) for f in self.func_dJay]
        self.dJaz = [f(self.q_list, self.dq_list) for f in self.func_dJaz]
        self.dJo = [f(self.q_list, self.dq_list) for f in self.func_dJo]
        return
    
    def update_maps_q_to_x(self,):
        """q -> xのpsi, J, dJを更新"""
        for i in range(7):
            name = 'q_to_x' + str(i)
            self.maps[name] = [
                self.T[2+i][0:3, 3:4],
                self.Jo[2+i],
                self.dJo[2+i],
            ]
        return
    
    def update_maps_x_GL_to_dis_x_to_goal(self, g, dg):
        """x_GL -> dis_x_to_goalのpsi, J, dJを更新"""
        
        z = g - self.T[-1][0:3, 3:4]
        
        dis = np.linalg.norm(z)  # GLからgoalまでの距離
        
        psi = dis
        J = -1/dis * z.T
        
        dz = dg - J @ self.q
        dJ = 1/(dis**2) * z.T - 1/dis * (dg - dz).T
        
        self.maps['x_GL_to_dis_x_to_goal'] = [psi, J, dJ]
        return
    
    def update_x_to_dis_x_to_obs(self, o, do):
        """x -> dis_x_to_obsのpsi, J, dJを更新"""
        
    
    
    def update_all_maps(self, q, dq):
        """"psi, J, dJを全て更新"""
        self.q = q
        self.dq = dq
        self.q_list = np.ravel(q).tolist()
        self.dq_list = np.ravel(dq).tolist()
        
        # for id in range(1):
        #     self.maps[id] = ()
        
        return
    
    
    def spesify_map(self, name):
        """タスク写像，そのヤコビ行列，ヤコビ行列の微分を指定"""
        psi, J, dJ = self.maps[name]
        return psi, J, dJ





if __name__ == '__main__':
    
    q = np.array([[1, 2, 3, 4, 5, 6, 7,]]).T
    dq = np.array([[1, 2, 3, 4, 5, 6, 7,]]).T
    
    hoge = Maps()
    hoge.update_all_maps(q, dq)
    hoge.update_all_baxter()
    hoge.update_maps_q_to_x()
    
