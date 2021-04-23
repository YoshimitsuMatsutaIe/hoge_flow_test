"""作成中"""
import rmp

class RMPalgebla(rmp.OriginalRMP):
    """RMP代数？"""
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.goal = kwargs.pop('goal')
        self.obs = kwargs.pop('obs')
        self.jl = kwargs.pop('jl')
    
    
    def push(self):
        
        return None
    
    def pull(self):
        return None
    
    def resolve(self):
        return None





if __name__ == "__main__":
    pass