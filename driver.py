import cv2
import time
import numpy as np
import gudhi as gd
import PyDMTGraph

import matplotlib as mpl
import matplotlib.pyplot as plt

def blur(im, n_itr=7):
    im_blur = np.copy(im)
    for i in range(n_itr):
        im_blur = cv2.GaussianBlur(im_blur, (3,3), 0)
    return im_blur

def sharpen(im, n_itr=3):

    kernel = 1/34 * np.array([[ 0,  0, -1,  0,  0],
                              [ 0, -1, -2, -1,  0],
                              [-1, -2, 17, -2, -1],
                              [ 0, -1, -2, -1,  0],
                              [ 0,  0, -1,  0,  0]])
    im_sharpen = np.copy(im)
    for i in range(n_itr):
        im_sharpen += np.abs(cv2.filter2D(im_sharpen, -1, kernel))
        # im_sharpen += cv2.filter2D(im_sharpen, -1, kernel)
    return im_sharpen

def sobel_sharpen(im, n_itr=3):

    vertical = 1/8 * np.array([[ 1, 0, -1],
                               [ 2, 0, -2],
                               [ 1, 0, -1]])

    horizontal = vertical.T

    im_sharpen = np.copy(im)
    for i in range(n_itr):
        im_sharpen += np.abs(cv2.filter2D(im_sharpen, -1, vertical)) + \
                      np.abs(cv2.filter2D(im_sharpen, -1, horizontal))
    return im_sharpen




def circle(n=500):
    n = 500
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dist = np.sqrt((i - n/2)**2 + (j - n/2)**2)
            A[i,j] = np.abs(n/4 - dist)
    return A

def plotGraph(ax, vertices, edges):
    lines = []
    for e in edges:
        v1 = vertices[e[0]]
        v2 = vertices[e[1]]
        # matplotlib expects vertex coordinates in column-major format
        # so we reverse the coordinates here
        c1 = (v1[1], v1[0])
        c2 = (v2[1], v2[0])
        lines.append([c1, c2])

    lc = mpl.collections.LineCollection(lines, colors="red")
    ax.add_collection(lc)


# img = circle()
img = cv2.imread("./GraphData/Thick7.png", cv2.IMREAD_GRAYSCALE).astype(np.double)
img = -img + 255
img = cv2.resize(img, (500, 500))
# img = (img > 240).astype(int) * 255.0
img = blur(img)
img = sharpen(img, n_itr=4)
img = blur(img)

start_time = time.time()

G = PyDMTGraph.DMTGraph(img)
# pairs = G.persistentIntervalsInDimension(0)
# print(pairs)

# TODO: compute homology
vertices, edges = G.computeGraph(0.1, 1)


end_time = time.time()

print(end_time - start_time)



fig, ax = plt.subplots()
ax.imshow(img)
plotGraph(ax, vertices, edges)
# plt.imshow(img)
plt.show()
