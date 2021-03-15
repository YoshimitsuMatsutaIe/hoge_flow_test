"""３ｄグラフ（）"""

"""3dグラフ作って回してるだけ
https://stackoverflow.com/questions/10374930/matplotlib-annotating-a-3d-scatter-plot
"""
# import numpy as np

# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt
# from numpy.random import rand

# from matplotlib import animation

# m = rand(3,3) # m is an array of (x,y,z) coordinate triplets

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# for i in range(len(m)): # plot each point + it's index as text above
#   x = m[i,0]
#   y = m[i,1]
#   z = m[i,2]
#   label = i
#   ax.scatter(x, y, z, color='b')
#   ax.text(x, y, z, '%s' % (label), size=20, zorder=1, color='k')

# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

# def animate(frame):
#   ax.view_init(30, frame/4)
#   plt.pause(.001)
#   return fig

# anim = animation.FuncAnimation(fig, animate, frames=200, interval=50)
# plt.show()



"""テキストが動く
https://sabopy.com/py/matplotlib-animation-48/
"""
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
import numpy as np

# Attaching 3D axis to the figure
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
# Lines to plot in 3D
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.plot(x, y, z, label='parametric curve')
ax.legend()

texts = []

def update(num):
    if len(texts) > 0:
        texts.pop().remove()
    text_, = [ax.text(x[num], y[num], z[num], "sabopy.com", color='red')]
    texts.append(text_)

ani = animation.FuncAnimation(fig, update, 100,interval=100)

#ani.save('text_anim.mp4', writer="ffmpeg",dpi=100)
plt.show()