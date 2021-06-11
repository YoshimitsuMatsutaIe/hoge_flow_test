"""rmp_mainとほぼ同じ"""


import numpy as np
import math
from math import pi, cos, sin, tan

import rmp_tree

def soft_normal(v, alpha):
    """ソフト正規化関数"""
    v_norm = np.linalg.norm(v)
    softmax = v_norm + 1 / alpha * np.log(1 + np.exp(-2 * alpha * v_norm))
    return v / softmax

def metric_stretch(v, alpha):
    """空間を一方向に伸ばす計量"""
    xi = soft_normal(v, alpha)
    return xi @ xi.T

def basic_metric_H(f, alpha, beta):
    """基本の計量"""
    f_norm = np.linalg.norm(f)
    f_softmax = f_norm + 1 / alpha * np.log(1 + np.exp(-2 * alpha * f_norm))
    s = f / f_softmax
    return beta * s @ s.T + (1 - beta) * np.eye(3)

def sigma_L(q, q_min, q_max):
    """アフィン変換andシグモイド変換写像"""
    return (q_max - q_min) * (1 / (1 + np.exp(-q))) + q_min

def D_sigma(q, q_min, q_max):
    """ジョイント制限に関する対角ヤコビ行列"""
    diags = (q_max - q_min) * (np.exp(-q) / (1 + np.exp(-q)) ** 2)
    #print(np.diag(diags.ravel()))
    return np.diag(diags.ravel())



class TargetAttractor(rmp_tree.RMPLeafBase):
    """アトラクター"""
    
    def __init__(
        self, name, parent, parent_param,
        attract_max_speed, attract_gain, attract_a_damp_r,
        attract_sigma_W, attract_sigma_H, attract_A_damp_r,
        robot_model,
    ):
        """"""
        
        self.attract_max_speed = attract_max_speed
        self.attract_gain = attract_gain
        self.attract_damp = attract_gain / attract_max_speed
        self.attract_a_damp_r = attract_a_damp_r
        self.attract_sigma_W = attract_sigma_W
        self.attract_sigma_H = attract_sigma_H
        self.attract_A_damp_r = self.attract_A_damp_r
        
        
        super().__init__(
            name = name,
            parent = parent,
            psi = robot_model.glipper_o,
            J = robot_model.J_glipper_o,
            dJ = robot_model.dJ_glipper_o,
            rmp = self.rmp,
        )
    
    
    def rmp(self, x, dx):
        
        a = self.attract_gain * soft_normal(x, self.attract_a_damp_r) - self.attract_damp * dx
        
        dis = np.linalg.norm(x)
        weight = math.exp(-dis / self.attract_sigma_W)
        beta_attract = 1 - math.exp(-1 / (2 * (dis / self.attract_sigma_H) ** 2))
        M = weight * basic_metric_H(a, self.attract_A_damp_r, beta_attract)
        
        f = M @ a
        
        return f, M