import matplotlib as mpl
import matplotlib.pyplot as plt

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

    lc = mpl.collections.LineCollection(lines, colors="red" )
    ax.add_collection(lc)
