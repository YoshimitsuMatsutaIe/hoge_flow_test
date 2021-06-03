"""メインのやつ
.baxter
"""

import numpy as np
import scipy.integrate as integrate
from math import cos, sin, tan, pi
#import itertools
import csv
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from mpl_toolkits.mplot3d import Axes3D 
#import matplotlib.patches as patches
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname=r'C:\WINDOWS\Fonts\Arial.ttf', size=14)
import time

# from baxter_utils_2 import *
# from baxter_utils import *
import baxter_utils_3
import rmp



class Simu:
    
    # 各部長さ
    L = 278e-3,
    h = 64e-3,
    H = 1104e-3,
    L0 = 270.35e-3,
    L1 = 69e-3,
    L2 = 364.35e-3,
    L3 = 69e-3,
    L4 = 374.29e-3,
    L5 = 10e-3,
    L6 = 368.3e-3
    
    # ジョイント制限
    q1_min, q1_max = -141, 51
    q2_min, q2_max = -123, 60
    q3_min, q3_max = -173, 173
    q4_min, q4_max = -3, 150
    q5_min, q5_max = -175, 175
    q6_min, q6_max = -90, 120
    q7_min, q7_max = -175, 175
    q_min = np.array([[q1_min, q2_min, q3_min, q4_min, q5_min, q6_min, q7_min]]).T * pi / 180
    q_max = np.array([[q1_max, q2_max, q3_max, q4_max, q5_max, q6_max, q7_max]]).T * pi / 180
    
    
    def __init__(self, goal, obs, q_initm, dq_init,):
        self.goal = goal
        self.obs = obs
        
    