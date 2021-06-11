"""作成中"""

import numpy as np


class RMPAlgebra:
    
    @classmethod
    def pushforward(cls, x, dx, psi, J):
        """1個だけpush"""
        y = psi(x)
        dy = J(x) @ dx
        return y, dy
    
    @classmethod
    def pullback(cls, x, dx, J, dJ, f, M):
        """1個だけpull"""
        fi = J(x).T @ (f - M @ dJ(x, dx) @ dx)
        Mi = J(x).T @ M @ J(x)
        return fi, Mi
    
    @classmethod
    def resolve(cls, f, M):
        a = np.linalg.pinv(M) @ f
        return a


class RMPNode:
    """RMPtreeのnode（rootとleaf以外）"""
    
    def __init__(self, name, parent, psi, J, dJ,):
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
        """子供を追加"""
        self.children.append(child)
        return
    
    
    def pushforward(self,):
        """parentノード ⇒ childノード
        
        再帰関数？
        """
        
        if self.psi is not None and self.J is not None:
            self.x, self.dx = RMPAlgebra.pushforward(self.dx, self.dx, self.psi, self.J)
        
        [child.pushforward() for child in self.children]
    
    
    def pullback(self,):
        """childノード -> parentノード"""
        
        [child.pullback for child in self.children]
        
        f = np.zeros_like(self.x)
        M = np.zeros((max(self.x.shape), max(self.x.shape)))
        
        for child in self.children:
            if child.f is not None and child.M is not None:
                fi, Mi = RMPAlgebra.pullback(
                    self.x, self.dx, child.J, child.dJ, child.f, child.M,
                    )
                f += fi
                M += Mi
        
        self.f = f
        self.M = M



class RMPRoot(RMPNode):
    """rootノード"""
    
    def __init__(self, name):
        RMPNode.__init__(
            self, name=name, parent=None, psi=None, J=None, dJ=None
        )
        return
    
    
    def set_root_state(self, x, dx,):
        self.x = x
        self.dx = dx
        return
    
    
    def pushforward(self):
        [child.pushforward for child in self.children]
    
    
    def resolve(self,):
        self.a = RMPAlgebra.resolve(self.f, self.M)
        return self.a
    
    
    def solve(self, x, dx,):
        self.set_root_state(x, dx)
        self.pushforward()
        self.pullback()
        optiomal_ddq = self.resolve()
        
        return optiomal_ddq



class RMPLeafBase(RMPNode):
    """leafノードのベース"""
    
    def __init__(self, name, parent, parent_param, psi, J, dJ, rmp):
        RMPNode.__init__(
            self, name=name, parent=parent, psi=psi,J=J, dJ=dJ
        )
        self.rmp = rmp
        self.parent_param = parent_param
    
    
    def eval_leaf(self,):
        self.f, self.M = self.rmp(self.x, self.dx)
        return
    
    
    def pullback(self,):
        self.eval_leaf()
        return
    
    
    def update(self):
        self.pushforward()
        return


if __name__ == "__main__":
    pass