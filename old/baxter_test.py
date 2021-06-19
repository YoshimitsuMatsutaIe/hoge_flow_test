"""baxter関連をテストする"""

import numpy as np
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



BaxterKinema = baxter_utils_3.BaxterKinematics3(
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
)
