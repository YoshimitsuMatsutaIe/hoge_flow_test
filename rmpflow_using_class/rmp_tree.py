"""作成中"""

import numpy as np
#import types

class RMPAlgebra:
    
    @classmethod
    def pushforward(cls, x, dx, psi, J):
        """1個だけpush"""
        y = psi
        dy = J @ dx
        return y, dy
    
    @classmethod
    def pullback(cls, x, dx, J, dJ, f, M):
        """1個だけpull"""
        fi = J.T @ (f - M @ dJ @ dx)
        Mi = J.T @ M @ J
        return fi, Mi
    
    @classmethod
    def resolve(cls, f, M):
        a = np.linalg.pinv(M) @ f
        return a


class RMPNode:
    """RMPtreeのnode（rootとleaf以外）"""
    
    def __init__(self, name, parent, psi=None, J=None, dJ=None, map_set=None):
        self.name = name
        self.parent = parent
        self.children = []
        
        # nodeをparnetに接続＿
        if self.parent:
            self.parent.add_child(self)
        
        if map_set:
            self.psi = map_set[0]
            self.J = map_set[1]
            self.dJ = map_set[2]
        else:
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
        
        if self.psi is not None:
            if self.J is not None:
                self.x, self.dx = RMPAlgebra.pushforward(self.dx, self.dx, self.psi, self.J)
        
        [child.pushforward() for child in self.children]
    
    
    def pullback(self,):
        """childノード -> parentノード"""
        
        [child.pullback for child in self.children]
        
        f = np.zeros_like(self.x)
        M = np.zeros((max(self.x.shape), max(self.x.shape)))
        
        for child in self.children:
            if child.f is not None:
                if child.M is not None:
                    fi, Mi = RMPAlgebra.pullback(
                        self.x, self.dx, child.J, child.dJ, child.f, child.M,
                        )
                    f += fi
                    M += Mi
        
        self.f = f
        self.M = M
    
    
    def remove_children(self,):
        """childrenを初期化
        ・障害物点の数が変動する場合に使用
        """
        
        self.children = []
        return



class RMPRoot(RMPNode):
    """rootノード"""
    
    def __init__(self, name):
        RMPNode.__init__(
            self, name=name, parent=None, psi=None, J=None, dJ=None, map_set=None,
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
    
    def __init__(self, name, parent, parent_param, rmp, psi=None, J=None, dJ=None, map_set=None,):
        RMPNode.__init__(
            self, name=name, parent=parent, psi=psi, J=J, dJ=dJ, map_set=map_set,
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