import sys
import imageio
from matplotlib.pyplot import figure, draw, pause, gca
import matplotlib.pyplot as plt
import numpy as np


def compareGraphs(u, v, im2, output):

    ax = figure().gca()
    ax.set_axis_off()
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_locator(plt.NullLocator())
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
                        hspace=0, wspace=0)
    ax.imshow(im2, cmap='gray', origin='lower')
    # plt.scatter(POI[:,0,1],POI[:,0,0])
    for i in range(0, u.shape[0], 10):
        for j in range(0, v.shape[1], 10):
            ax.arrow(
                j,
                i,
                v[i, j] * 4,
                u[i, j] * 4,
                color='red',
                head_width=0.1,
                head_length=0.1,
            )

    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    draw()
    pause(0.01)
    plt.savefig(output, bbox_inches='tight',
                pad_inches=0, dpi=200)


u_tmp = []
v_tmp = []
f = open(sys.argv[1])
for line in f:
    t = line.split()
    t = [float(x) for x in t]
    u_tmp.append(t)
f = open(sys.argv[2])
for line in f:
    t = line.split()
    t = [float(x) for x in t]
    v_tmp.append(t)
u_tmp = np.array(u_tmp)
v_tmp = np.array(v_tmp)

im1 = imageio.imread(sys.argv[3], as_gray=True)
im2 = imageio.imread(sys.argv[4], as_gray=True)
compareGraphs(u_tmp, v_tmp, im2, sys.argv[5])
