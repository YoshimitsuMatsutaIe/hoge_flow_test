

class Tasu:
    def __init__(self, a):
        self.a = a
    
    def tasu(self, b):
        return self.a + b


class Hiku:
    def __init__(self, a):
        self.a = a
    
    def hiku(self, b):
        return self.a - b


class TasuHiku(Tasu, Hiku):
    def __init__(self):
        