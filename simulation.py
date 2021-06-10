"""シミュレーションを実行

scupy.integrateを使用
"""

import numpy as np
import scipy.integrate as integrate


class Simulation:
    """"""
    
    def __init__(
        self,
        rob_model, rmp, TIME_SPAN=10, TIME_INTERVAL=0.01,
    ):
        self.TIME_SPAN = TIME_SPAN
        self.TIME_INTERVAL = TIME_INTERVAL
        self.rob_model = rob_model
        self.rmp = rmp
    
    
    def run_simu(self, q0, goal, obs):
        """"""
        
        def diff_eq(t, x,):
            q = np.array([x[0:len(x)/2]]).T
            dq = np.array([x[len(x)/2:]]).T
            
            ddq = self.rmp.compute_optiomal_ddq(q, dq, goal, obs)
            
            return dq, ddq
        
        sol = integrate.solve_ivp(
            
        )