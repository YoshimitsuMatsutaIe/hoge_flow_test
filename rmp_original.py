"""rmp_mainとほぼ同じ"""


import numpy as np
import math
from math import pi, cos, sin, tan
import autograd.numpy as agnp
import autograd

import rmp_tree

def soft_normal(v, alpha):
    """ソフト正規化関数"""
    v_norm = np.linalg.norm(v)
    softmax = v_norm + 1 / alpha * math.log(1 + math.exp(-2 * alpha * v_norm))
    return v / softmax

def metric_stretch(v, alpha):
    """空間を一方向に伸ばす計量"""
    xi = soft_normal(v, alpha)
    return xi @ xi.T

def basic_metric_H(f, alpha, beta):
    """基本の計量"""
    f_norm = np.linalg.norm(f)
    f_softmax = f_norm + 1 / alpha * math.log(1 + math.exp(-2 * alpha * f_norm))
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



class TargetAttractorFromOriginal(rmp_tree.RMPLeafBase):
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


class CollisionAvoidanceFromOriginal(rmp_tree.RMPLeafBase):
    """障害物回避"""
    
    def __init__(
        self, name, parent, parent_param,
        obs_scale_rep, obs_scale_damp, obs_ratio, obs_rep_gain, obs_r,
        robot_model,
    ):
        self.obs_scale_rep = obs_scale_rep
        self.obs_scale_damp = obs_scale_damp
        self.obs_ratio = obs_ratio
        self.obs_rep_gain = obs_rep_gain
        self.obs_r = obs_r
        
        # なおしてほしい↓
        super().__init__(
            name = name,
            psi = robot_model.glipper_o,
            J = robot_model.J_glipper_o,
            dJ = robot_model.dJ_glipper_o,
            rmp = self.rmp,
        )
    
    
    def rmp(self, x, dx):
        damp_gain = self.obs_rep_gain * self.obs_ratio
        
        dis = np.linalg.norm(x)
        grad_dis = x / dis
        
        # 斥力項．障害物に対する位置ベースの反発力？
        alpha_rep = self.obs_rep_gain * math.exp(-dis / self.obs_scale_rep)  # 斥力の活性化関数
        a_rep = alpha_rep * grad_dis  # 斥力項
        
        # ダンピング項．障害物に向かう速度にペナルティを課す？
        P_obs = max(0, -(dx).T @ grad_dis) * grad_dis @ grad_dis.T @ dx  # 零空間射影演算子
        alpha_damp = damp_gain / (dis / self.obs_scale_damp + 1e-7)  # ダンピングの活性化関数
        a_damp = alpha_damp * P_obs  # ダンピング項
        #print("a_obs = ", a_rep + a_damp)
        a =  a_rep + a_damp
        
        
        weight_obs = (dis / self.obs_r) ** 2 - 2 * dis / self.obs_r + 1
        M = weight_obs * np.eye(3, dtype = np.float32)
        
        f = M @ a
        
        return f, M



class TargetAttracttorFromGDS(rmp_tree.RMPLeafBase):
    """GDSからのやつ"""
    
    def __init__(
        self, name, parent, parent_param,
        attract_max_speed, attract_gain, attract_alpha_f,
        attract_sigma_alpha, attract_sigma_gamma,
        attract_w_upper, attract_w_lower,
        attract_alpha, attract_epsilon,
        robot_model,
    ):
        self.attract_max_speed = attract_max_speed
        self.attract_gain = attract_gain
        self.attract_alpha_f = attract_alpha_f
        self.attract_sigma_alpha = attract_sigma_alpha
        self.attract_sigma_gamma = attract_sigma_gamma
        self.attract_w_upper = attract_w_upper
        self.attract_w_lower = attract_w_lower
        self.attract_alpha = attract_alpha
        self.attract_epsilon = attract_epsilon
        
        super().__init__(
            name = name,
            parent = parent,
            J = None,
            dJ = None,
            rmp = self.rmp,
        )
    
    
    def rmp(self, x, dx):
        
        norm_x = np.linalg.norm(x)
        hat_x = x / norm_x
        
        # 慣性行列を計算
        alpha = math.exp(-norm_x**2 / (2*self.attract_sigma_alpha**2))
        gamma = math.exp(-norm_x**2 / (2*self.attract_sigma_gamma**2))
        w = gamma * self.attract_w_upper + (1 - gamma) * self.attract_w_lower
        
        def s_alpha(alpha, hat_x):
            return (1 - math.exp(-2*alpha*hat_x)) / (1 + math.exp(-2*alpha*hat_x))
        
        nabla_x_potential = s_alpha(self.attract_alpha, hat_x) * hat_x
        
        M = w * ((1 - alpha) * (nabla_x_potential @ nabla_x_potential.T) + \
            (alpha + self.attract_epsilon) * np.eye(3))
        
        
        # 曲率項xiを計算
        def M_stretch(x):
            norm_x = agnp.linalg.norm(x)
            hat_x = x / norm_x
            
            alpha = agnp.exp(-norm_x**2 / (2*self.attract_sigma_alpha**2))
            gamma = agnp.exp(-norm_x**2 / (2*self.attract_sigma_gamma**2))
            w = gamma * self.attract_w_upper + (1 - gamma) * self.attract_w_lower
            
            def s_alpha(alpha, hat_x):
                return (1 - agnp.exp(-2*alpha*hat_x)) / (1 + agnp.exp(-2*alpha*hat_x))
            
            nabla_x_potential = s_alpha(self.attract_alpha, hat_x) * hat_x
            
            M = w * ((1 - alpha) * (nabla_x_potential @ nabla_x_potential.T) + \
                (alpha + self.attract_epsilon) * agnp.eye(3))
            
            return M
        
        def partial_mi_s(x):
            
            x = agnp([[x[0,0], x[1,0], x[2, 0]]])
            
            def m1(x):
                return M_stretch(x)[:, 0]
            
            def m2(x):
                return M_stretch(x)[:, 1]
            
            def m3(x):
                return M_stretch(x)[:, 2]
            
            partial_x_m1 = autograd.jacobian(m1, x)
            partial_x_m2 = autograd.jacobian(m2, x)
            partial_x_m3 = autograd.jacobian(m3, x)
            
            return [partial_x_m1, partial_x_m2, partial_x_m3]
        
        partial_ms = partial_mi_s(x)
        xi_1_ = [par @ dx for par in partial_ms]
        xi_1 = np.concatenate(xi_1_, axis = 1) @ dx
        
        def calc_xi_2(x, dx):
            def dxT_M_dx(x):
                z = dx.T @ M_stretch(x) @ dx
                return np.ravel(z)
            
            x = agnp([[x[0,0], x[1,0], x[2, 0]]])
            xi_2 = 1/2 * autograd.grad(dxT_M_dx)
            
            return xi_2
        
        
        xi_2 = calc_xi_2(x, dx)
        
        xi_M = xi_1 - xi_2
        
        
        # 力を計算
        gamma_p = self.attract_gain
        gamma_d = self.attract_gain / self.attract_max_speed
        
        f_ = -gamma_p * soft_normal(x, dx) - gamma_d * dx
        carv_ = -np.linalg.inv(M) @ xi_M
        f = f_ + carv_
        
        return f, M







if __name__ == '__main__':
    pass