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
    softmax = v_norm + 1 / alpha * np.log(1 + np.exp(-2 * alpha * v_norm))
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



class TargetAtracttorFromGDS(rmp_tree.RMPLeafBase):
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
        
        # 暫定的
        super().__init__(
            name = name,
            parent = parent,
            parent_param=None,
            psi = None,
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
        
        def s_alpha(alpha, norm_x):
            return (1 - math.exp(-2*alpha*norm_x)) / (1 + math.exp(-2*alpha*norm_x))
        
        nabla_x_potential = s_alpha(self.attract_alpha, norm_x) * hat_x
        
        M = w * ((1 - alpha) * (nabla_x_potential @ nabla_x_potential.T) + \
            (alpha + self.attract_epsilon) * np.eye(3))
        
        
        # 曲率項xiを計算
        def M_stretch(x):
            """autogradで使う用"""
            x = x.reshape(3, 1)
            norm_x = agnp.linalg.norm(x)
            hat_x = x / norm_x
            alpha = agnp.exp(-norm_x**2 / (2*self.attract_sigma_alpha**2))
            gamma = agnp.exp(-norm_x**2 / (2*self.attract_sigma_gamma**2))
            w = gamma * self.attract_w_upper + (1 - gamma) * self.attract_w_lower
            
            def s_alpha(alpha, norm_x):
                return (1 - agnp.exp(-2*alpha*norm_x)) / (1 + agnp.exp(-2*alpha*norm_x))
            
            nabla_x_potential = s_alpha(self.attract_alpha, norm_x) * hat_x
            
            M = w * ((1 - alpha) * (nabla_x_potential @ nabla_x_potential.T) + \
                (alpha + self.attract_epsilon) * np.eye(3))
            return M
        
        def partial_mi_s(x):
            """autogradで面倒な計算をやる"""
            def m1(x):
                return M_stretch(x)[:, 0]
            
            def m2(x):
                return M_stretch(x)[:, 1]
            
            def m3(x):
                return M_stretch(x)[:, 2]
            
            partial_x_m1 = autograd.jacobian(m1)(x)
            partial_x_m2 = autograd.jacobian(m2)(x)
            partial_x_m3 = autograd.jacobian(m3)(x)
            
            return [partial_x_m1, partial_x_m2, partial_x_m3]
        
        partial_ms = partial_mi_s(np.ravel(x))
        xi_1_ = [par @ dx for par in partial_ms]
        xi_1 = np.concatenate(xi_1_, axis = 1) @ dx
        #print('xi_1=\n', xi_1)
        
        def calc_xi_2(x, dx):
            def dxT_M_dx(x):
                z = dx.T @ M_stretch(x) @ dx
                return np.ravel(z)
            
            x = np.ravel(x)
            xi_2 = 1/2 * autograd.grad(dxT_M_dx)(x)
            
            return xi_2.reshape(3, 1)
        
        
        xi_2 = calc_xi_2(x, dx)
        #print('xi_2=\n', xi_2)
        
        xi_M = xi_1 - xi_2
        #print('xi_M=\n', xi_M)
        
        # 力を計算
        gamma_p = self.attract_gain
        gamma_d = self.attract_gain / self.attract_max_speed
        
        f_ = M @ (-gamma_p * soft_normal(x, dx) - gamma_d * dx)
        #print(f_)
        carv_ = xi_M
        #print(carv_)
        f = f_ + carv_
        
        return f, M



class CollisionAvoidanceFromGDS(rmp_tree.RMPLeafBase):
    """GDSを使った衝突回避RMP"""
    
    def __init__(
        self, name, parent, parent_param,
        obs_rw, obs_sigma, obs_alpha,
        
    ):
        
        self.obs_rw = obs_rw
        self.obs_sigma = obs_sigma
        self.obs_alpha = obs_alpha
        
        super().__init__(
            name = name,
            parent = parent,
            parent_param=None,
            psi = None,
            J = None,
            dJ = None,
            rmp = self.rmp,
        )
    
    
    def rmp(self, x, dx):
        """1d"""
        
        def w(s):
            z = max(self.obs_rw - s, 0)
            return z**2 / s
        
        def u(ds):
            if ds < 0:
                return 1 - np.exp(-ds**2 / (2 * self.obs_sigma**2))
            else:
                return 0
        
        def dw(s):
            if s < self.obs_rw:
                return -(self.obs_rw - s)**2 / s**2 - 2*(self.obs_rw - s) / s
            else:
                return 0
        
        def du(ds):
            if ds < 0:
                return -ds**3 / (2*self.obs_sigma**4) * np.exp(-ds**2 / (2*self.obs_sigma**2))
            else:
                return 0
        
        def grad_potential(s):
            return self.obs_alpha * w(s) * dw(s)
        
        # 慣性行列
        M = w(x) * (u(dx) + 1/2 * dx * du(dx))
        
        # 力
        xi = 1/2 * u(dx) * dw(x) * dx**2
        #grad_tilda_potential = 1/w(x) * grad_potential(x)  # fのときw消えるから無意味
        f = -grad_potential(x) - xi
        
        return f, M


class JointLimitAvoidanceFromGDS(rmp_tree.RMPLeafBase):
    """ジョイント制限"""
    
    def __init__(
        self, name, parent, parent_param,
        jl_c, jl_lambda, jl_upper, jl_lower, jl_eta_p, jl_eta_d,
    ):
        
        self.jl_c = jl_c
        self.jl_lambda = jl_lambda
        self.jl_upper = jl_upper
        self.jl_lower = jl_lower
        self.jl_eta_p = jl_eta_p
        self.jl_eta_d = jl_eta_d
        
        super().__init__(
            name = name,
            parent = parent,
            parent_param=None,
            psi = None,
            J = None,
            dJ = None,
            rmp = self.rmp,
        )
    
    
    def rmp(self, q, dq):
        
        # 速度依存計量を計算
        def sigma_i(u):
            return 1 / (1 + np.exp(-u))
        
        def d_ii(i, u):
            return (self.jl_upper[i] - self.jl_lower[i]) * (1 + np.exp(-u))**(-2) * np.exp(-u)
        
        def alpha_i(du):
            return 1 / (1 + np.exp(-self.jl_c * du))
        
        def tilda_dii(i, u, du):
            dii = d_ii(i, u)
            sigma = sigma_i(u)
            alpha = alpha_i(du)
            return sigma * (alpha*dii + (1-alpha)) + (1-sigma) * ((1-alpha)*dii + alpha)
        
        def tilda_D_sigma(q, dq):
            Dii = [tilda_dii(i, q[i, 0], dq[i, 0]) for i in range(q.shape[0])]
            return np.diag(Dii)
        
        inv_tilda_D_sigma = np.linalg.inv(tilda_D_sigma(q, dq))
        A = self.jl_lambda * np.linalg.matrix_power(inv_tilda_D_sigma, 2)
        
        
        # 力を計算
        def dsigma_idq_i(u):
            return (1 + np.exp(-u))**(-2) * np.exp(-u)
        
        def dd_iidq_i(i, u):
            return (self.jl_upper[i] - self.jl_lower[i]) * \
                (2 * (1 + np.exp(-u))**(-3) * np.exp(-2*u) + \
                    (1 + np.exp(-u))**(-2) * -np.exp(-u))
        
        def dAiidqi(i, u, du):
            dii = d_ii(i, u)
            sigma = sigma_i(u)
            alpha = alpha_i(du)
            dsigma = dsigma_idq_i(u)
            ddii = dd_iidq_i(i, u)
            return dsigma * (alpha*dii + (1-alpha)) + \
                sigma * (alpha*ddii) + \
                    -dsigma * ((1-alpha)*dii + alpha) + \
                        (1-sigma) * ((1-alpha)*ddii)
        
        xi_A_ii = [1/2 * dAiidqi(i, q[i,0], dq[i,0]) * dq[i,0]**2 for i in range(q.shape[0])]
        xi_A = np.array([xi_A_ii]).T
        ddq = A @ (self.jl_eta_p * (-q) - self.jl_eta_d * dq) - xi_A
        
        return ddq, A


# テスト
def test_TaregetAtracttorFromGDS():
    """test"""
    hoge = TargetAtracttorFromGDS(
        name='hoge', parent = None, parent_param=None,
        attract_max_speed = 0.1, 
        attract_gain = 1,
        attract_alpha_f = 0.3,
        attract_sigma_alpha = 1,
        attract_sigma_gamma = 1,
        attract_w_upper = 1,
        attract_w_lower = 1,
        attract_alpha = 0.1,
        attract_epsilon = 1e-5,
        robot_model = None,
    )
    
    x = np.array([[1., 1., 1.]]).T
    dx = np.array([[1., 13., 1.]]).T
    
    f, M = hoge.rmp(x, dx)
    print('f=\n', f)
    print('M=\n', M)
    return


def test_CollisionAvoidanceFromGDS():
    hoge = CollisionAvoidanceFromGDS(
        name='hoge', parent = None, parent_param=None,
        obs_rw = 1,
        obs_sigma=1,
        obs_alpha=1,
    )
    
    x = 0.1
    dx = -0.5
    
    f, M = hoge.rmp(x, dx)
    print('f=\n', f)
    print('M=\n', M)
    return


def test_JointLimitAvoidanceFromGDS():
    jl_upper = [pi] * 7
    jl_lower = [-pi] * 7
    hoge = JointLimitAvoidanceFromGDS(
        name='hoge', parent = None, parent_param=None,
        jl_c = 1, jl_lambda = 1, jl_upper = jl_upper, jl_lower= jl_lower,
        jl_eta_p=1, jl_eta_d=1,
    )
    
    q = np.array([[0, 0, 0, 0, 0, 0, 0]]).T
    dq = np.array([[0, 0, 0, 0, 0, 0, 0]]).T
    
    f, M = hoge.rmp(q, dq)
    print('f=\n', f)
    print('M=\n', M)
    
    return



if __name__ == '__main__':
    #test_CollisionAvoidanceFromGDS()
    
    test_JointLimitAvoidanceFromGDS()
