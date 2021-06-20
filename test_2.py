import numpy as np
import types

def a():
    return 2

print(type(a))

b = 1

if isinstance(a, types.FunctionType):
    print('関数')
else:
    pass