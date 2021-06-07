"""作成中"""

import numpy as np

import rmp



class RMPAlgebla:
    """RMP代数"""
    
    def __init__(self):
        return
    
    
    def push(self, psi, J):
        """pushforword演算"""
        
        if self.psi is not None and self.J is not None:
            self.x = self.psi(self.parent.x)
            self.dx = self.J(self.parent.x) @ self.parent.dx
        
        [child.push() for child in self.children]
    
    
    def pull(self,):
        """pullback演算"""
        
        [child.pull() for child in self.children]
        
        f = np.zeros_like(self.x)
        M = np.zeros((max(self.x.shape), max(self.x.shape)))
        
        for child in self.children:
            J_child = child.J(self.x)
            dJ_child = child.dJ(self.x, self.dx)
            
            if child.f is not None and child.M is not None:
                f += J_child.T @ \
                    (child.f - child.M @ dJ_child @ self.dx)
                M += J_child.T @ child.M @ J_child
        
        self.f = f
        self.M = M
    
    
    def resolve(self):
        return None



class RMPNode(RMPAlgebla):
    """RMPtreeのnode"""
    
    def __init__(self, name, parent, psi, J, dJ,):
        
        RMPAlgebla.__init__(self,)
        
        self.name = name
        self.parent = parent
        self.children = []
        self.psi = psi
        self.J = J
        self.dJ = dJ
        
        self.x = None
        self.dx = None
        
        self.f = None
        self.a = None
        self.M = None
    
    
    def add_child(self, child):
        self.children.append(child)


class RMPTree:
    
    def __init__(self,):
        pass








if __name__ == "__main__":
    pass