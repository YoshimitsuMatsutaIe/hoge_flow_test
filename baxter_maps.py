"""baxterのタスク写像"""


import numpy as np

import baxter_oldold


class Maps(baxter_oldold.BaxterFuncs):
    """baxterの写像"""
    
    def __init__(self,):
        
        super().__init__()
        self.init_r_bars()
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
    
    
    def maps_q_to_x(self,):
        self.maps['q_to_x0'] = ()
    
    
    def update_all_maps(self, q, dq):
        """"psi, J, dJを全て更新"""
        
        self.q = np.ravel(q).tolist()
        self.dq = np.ravel(dq).tolist()
        
        for id in range(1):
            self.maps[id] = ()
        
        return
    
    
    def spesify_map(self, name):
        """タスク写像，そのヤコビ行列，ヤコビ行列の微分を指定"""
        psi, J, dJ = self.maps[name]
        return psi, J, dJ





if __name__ == '__main__':
    
    hoge = Maps()