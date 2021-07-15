import numpy as np
#import types


def pushforward(x, dx, psi, J):
    """1個だけpush"""
    y = psi
    dy = J @ dx
    return y, dy
    

def pullback(x, dx, J, dJ, f, M):
    """1個だけpull"""
    fi = J.T @ (f - M @ dJ @ dx)
    Mi = J.T @ M @ J
    return fi, Mi
    

def resolve(f, M):
    a = np.linalg.pinv(M) @ f
    return a



class RMPNode:
    """ノード"""
    
    def __init__(self, name, parent,):
        self.name = name
        self.parent = parent
        self.children = []
        
        # 自分をparent-nodeに接続
        if self.parent:
            self.parent.add_child(self)
        
        self.x = None
        self.dx = None
        self.f = None
        self.a = None
        self.M = None
        
        return
    
    
    def add_child(self, child):
        self.children.append(child)
        return
    
    
    def pushforward(self,):
        
        
        
    
    