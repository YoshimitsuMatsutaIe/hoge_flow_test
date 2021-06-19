"""メインのやつ
.baxter
"""

import numpy as np
from scipy.integrate import solve_ivp
from math import cos, sin, tan, pi
#import itertools
import csv
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from mpl_toolkits.mplot3d import Axes3D 
#import matplotlib.patches as patches
import time

# 自作
import rmp_tree
import rmp_leaf


def simu():

    q1_init = 0
    q2_init = -31
    q3_init = 0
    q4_init = 43
    q5_init = 0
    q6_init = 72
    q7_init = 0
    q = np.array([[q1_init, q2_init, q3_init, q4_init, q5_init, q6_init, q7_init]]).T * pi / 180 # ジョイント角度ベクトル
    dq = np.array([[0, 0, 0, 0, 0, 0, 0]]).T  # ジョイント角速度ベクトル
    
    state_init = np.concatenate((q, dq), axis=None)
    #print(state_init)
    
    
    ### tree構築 ###
    # root
    root = rmp_tree.RMPRoot(name='root')
    
    # frame
    frames = []
    for i in range(7):
        name = 'T_' + str(i)
        frames.append(
            rmp_tree.RMPNode(name, root)
        )
    
    
    
    
    def dynamics(t, state):
        state = state.reshape((2, 1))
        q = state[0]
        dq = state[1]
        
        ddq = root.solve(x=q, dx=dq)
        
        return np.concatenate((dq, ddq), axis = None)
    
    TIME_SPAN = 5.00
    TIME_INTERVAL = 0.01
    t = np.arange(0, TIME_SPAN, TIME_INTERVAL)
    
    ### シミュレーション実行 ###
    sol = solve_ivp(
        fun = dynamics,
        t_span = (0, TIME_SPAN),
        y0 = state_init,
        method = 'RK45',
        t_eval = t,
    )
    
    
    ### 図示 ###
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    plt.show()




if __name__ == '__main__':
    simu()