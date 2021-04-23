"""これが一番最新"""


import numpy as np
from math import cos, sin, tan, pi, sqrt



class BaxterKinematicsNew:
    """ローカル座標系の原点，ヤコビ行列等を計算
    ・
    """
    
    def __init__(self, **kwargs):
        """"""
        self.L = kwargs.pop('L')
        self.h = kwargs.pop('h')
        self.H = kwargs.pop('H')
        self.L0 = kwargs.pop('L0')
        self.L1 = kwargs.pop('L1')
        self.L2 = kwargs.pop('L2')
        self.L3 = kwargs.pop('L3')
        self.L4 = kwargs.pop('L4')
        self.L5 = kwargs.pop('L5')
        self.L6 = kwargs.pop('L6')
    
    def T_Wo_BL