"""シミュレーションを実行

scupy.integrateを使用
"""

from simu_main import Simu
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt 
import matplotlib.animation as anm

import baxter_utils_new
import rmp
import rmp_main


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
    
    
    def run_simu(self, q0, dq0, goal, obs):
        """"""
        
        def diff_eq(t, x,):
            q = np.array([x[0:len(x)/2]]).T
            dq = np.array([x[len(x)/2:]]).T
            
            ddq = self.rmp.compute_optiomal_ddq(q, dq, goal, obs)
            
            return dq, ddq
        
        sol = integrate.solve_ivp(
            fun = diff_eq,
            t_span = (0.0, self.TIME_SPAN),
            y0 = q0.append(dq0),
            method = 'RK45',
            t_eval = np.arange(0.0, self.TIME_SPAN, self.TIME_INTERVAL),
            #rtol = 1.e-12,
            #atol = 1.e-14,
        )
        
        return sol
    
    
    def draw(self, sol, goal, obs, anisave=False):
        """"""
        
        
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.grid(True)
        ax.set_xlabel('x[m]')
        ax.set_ylabel('y[m]')
        ax.set_zlabel('z[m]')
        
        ax.scatter(
            goal[:, 0], goal[:, 1], goal[:, 2],
        s = 100, label = 'goal point', marker = '*', color = '#ff7f00',
        alpha = 1, linewidths = 1.5, edgecolors = 'red',
        )
        ax.scatter(
        obs[:, 0], obs[:, 1], obs[:, 2],
        s = 100, label = 'obstacle point', marker = '+', color = 'k', 
        alpha = 1
        )
        
        



def main():
    
    model = baxter_utils_new.BaxterKinematicsNew()
    rmp = rmp_main.RMPTree()
    simu = Simulation()
    
    return

if __name__ == '__main__':
    main()