import numpy as np
import types



a = np.random.randint(1, 10, (3, 7))

print(a)

b = np.eye(7)
print(b)

print(a @ b)