import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
import imageio


X = np.arange(0, 26, 1)
Y = np.arange(0, 26, 1)
X, Y = np.meshgrid(X, Y)

fig = plt.figure()
ax = Axes3D(fig)
ax.view_init(25, 30)

for i in range(0, 60):
    Z = genfromtxt('csv/%s.csv' % (i * 10), delimiter=',')

    my_cmap = cm.coolwarm
    my_cmap.set_under('green', 0)
    ax.plot_surface(X, Y, Z, cmap=my_cmap, linewidth=0, vmin=0.2)
    plt.title('Анисимов (%s сек.)' % i)
    ax.set_zlabel('температура')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set(zlim=(0, 2000))
    filename = 'frames/' + str(i) + '.png'
    plt.savefig(filename, dpi=100)
    plt.cla()
    plt.gca()


with imageio.get_writer('graph_2.gif', mode='I') as writer:
    for i in range(0, 60):
        image = imageio.imread('frames/' + str(i) + '.png')
        os.remove('frames/' + str(i) + '.png')
        writer.append_data(image)
