
import numpy as np

a = np.array([[1, 2, 3, 4, 5]]).T

b = np.array([[2, 10, -5, 3, 0]]).T

print((a > b).any())
