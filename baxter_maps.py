"""baxterのタスク写像"""


import numpy as np
import time

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
        
        self.obs_num_old = None
        self.maps = {}


    def init_r_bars(self,):
        """制御点座標を追加"""

        r_bar_1 = [
            np.array([[0, self.L1/2, -self.L0/2, 1]]).T,
            np.array([[0, -self.L1/2, -self.L0/2, 1]]).T,
            np.array([[self.L1/2, 0, -self.L0/2, 1]]).T,
            np.array([[-self.L1/2, 0, -self.L0/2, 1]]).T,
        ]  # 1座標系からみた制御点位置

        r_bar_2 = [
            np.array([[0, 0, self.L3/2, 1]]).T,
            np.array([[0, 0, -self.L3/2, 1]]).T,
        ]

        r_bar_3 = [
            np.array([[0, self.L3/2, -self.L2*2/3, 1]]).T,
            np.array([[0, -self.L3/2, -self.L2*2/3, 1]]).T,
            np.array([[self.L3/2, 0, -self.L2*2/3, 1]]).T,
            np.array([[-self.L3/2, 0, -self.L2*2/3, 1]]).T,
            np.array([[0, self.L3/2, -self.L2*1/3, 1]]).T,
            np.array([[0, -self.L3/2, -self.L2*1/3, 1]]).T,
            np.array([[self.L3/2, 0, -self.L2*1/3, 1]]).T,
            np.array([[-self.L3/2, 0, -self.L2*1/3, 1]]).T,
        ]

        r_bar_4 = [
            np.array([[0, 0, self.L3/2, 1]]).T,
            np.array([[0, 0, -self.L3/2, 1]]).T,
        ]

        r_bar_5 = [
            np.array([[0, self.L5/2, -self.L4/2, 1]]).T,
            np.array([[0, -self.L5/2, -self.L4/2, 1]]).T,
            np.array([[self.L5/2, 0, -self.L4/2, 1]]).T,
            np.array([[-self.L5/2, 0, -self.L4/2, 1]]).T,
        ]

        r_bar_6 = [
            np.array([[0, 0, self.L5/2, 1]]).T,
            np.array([[0, 0, -self.L5/2, 1]]).T,
        ]

        r_bar_7 = [
            np.array([[0, self.L5/2, self.L6/2, 1]]).T,
            np.array([[0, -self.L5/2, self.L6/2, 1]]).T,
            np.array([[self.L5/2, 0, self.L6/2, 1]]).T,
            np.array([[-self.L5/2, 0, self.L6/2, 1]]).T,
        ]

        r_bar_GL = [
            np.array([[0, 0, 0, 1]]).T
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


    def update_maps_q_to_oGL(self,):
        """q -> oGLのpsi, J, dJを更新"""
        
        key = (-1)  # 0にしとく
        
        psi = self.T[-1][0:3, 3:4]
        J = self.Jo[-1]
        dJ = self.dJo[-1]

        self.maps[key] = [psi, J, dJ]
        return
    
    
    def update_q_to_x(self,):
        """q -> xのpsi, J, dJを更新"""
        
        for i, T in enumerate(self.T[2:]):  # 0, 1は固定
            
            Block_eye = np.concatenate([np.eye(3), np.zeros((3, 1))], axis = 1)
            Block_J = np.concatenate(
                [self.Jax[i], self.Jay[i], self.Jaz[i], self.Jo[i]],
                axis=1
            )  # ブロックヤコビ
            Block_dJ = np.concatenate(
                [self.dJax[i], self.dJay[i], self.dJaz[i], self.dJo[i]],
                axis=1
            )  # ブロック（微分）ヤコビ
            
            for j, r_bar in enumerate(self.r_bars[i]):
                
                key = (i, j)
                
                Block_r_bar = np.concatenate(
                    [
                        np.eye(7) * r_bar[0, 0],
                        np.eye(7) * r_bar[1, 0],
                        np.eye(7) * r_bar[2, 0],
                        np.eye(7),
                    ],
                    axis = 0,
                )
                
                psi = Block_eye @ T @ r_bar
                J = Block_J @ Block_r_bar
                dJ = Block_dJ @ Block_r_bar
                
                self.maps[key] = [psi, J, dJ]
        
        return
    
    
    
    def update_x_to_dis_x_to_obs(self, obs, dobs):
        """x -> dis_x_to_obsのpsi, J, dJを更新"""
        
        def maps_of_distance(c_point, dc_point, obj, dobj,):
            """c_point -> objectへの距離関数の汎用的な写像"""
            
            z = obj - c_point
            dis_z = np.linalg.norm(z)
            psi = dis_z
            J = -1/dis_z * z.T
            dz = dobj - J @ dc_point
            dJ = 1/(dis_z**2) * z.T - 1/dis_z * (dobj - dz).T
            
            return [psi, J, dJ]
        
        
        for i in range(7):
            for j in range(len(self.r_bars[i])):
                # 古いものを削除
                if self.obs_num_old is None:
                    pass
                else:
                    for k in range(self.obs_num_old):
                        key = (i, j, k)
                        del self.maps[key]
                
                for k, (o, do) in enumerate(zip(obs, dobs)):
                    key = (i, j, k)
                    
                    x, J_, _ = self.maps[(i, j)]
                    dx = J_ @ self.q
                    
                    self.maps[key] = maps_of_distance(x, dx, o, do)
        
        self.obs_num_old = len(obs)  # 1step前の障害物数の数
        return
    
    
    
    def update_all_maps(self, q, dq, obs, dobs):
        """"psi, J, dJを全て更新"""
        self.q = q
        self.dq = dq
        self.q_list = np.ravel(q).tolist()
        self.dq_list = np.ravel(dq).tolist()
        
        self.update_all_baxter()
        
        self.update_maps_q_to_oGL()
        self.update_q_to_x()
        self.update_x_to_dis_x_to_obs(obs, dobs)
        
        return


    def spesify_map(self, name):
        """タスク写像，そのヤコビ行列，ヤコビ行列の微分を指定"""
        psi, J, dJ = self.maps[name]
        return psi, J, dJ





if __name__ == '__main__':

    q = np.array([[1, 2, 3, 4, 5, 6, 7,]]).T
    dq = np.array([[1, 2, 3, 4, 5, 6, 7,]]).T
    
    o = [
        np.array([[1, 1, 1]]).T,
        np.array([[2, 5, 10]]).T,
    ]
    do = [
        np.array([[1, 1, 1]]).T,
        np.array([[2, 5, 10]]).T,
    ]
    
    start = time.time()
    hoge = Maps()
    for i in range(100):
        hoge.update_all_maps(q, dq, o, do)
    print('time=', time.time() - start)
    #print(hoge.maps)

