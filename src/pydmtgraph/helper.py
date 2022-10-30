import matplotlib as mpl
import matplotlib.pyplot as plt

def plotGraph(vertices, edges, ax):
    '''
        Plots the graph on the provided axis `ax`

        This method assumes the arrays `vertices` and
        `edges` have the same format as the arrays returned
        by the pydmtgraph.DMTGraph.computeGraph method.

        Args:
            vertices: a list of vertex positions,
            edges: a list of edge endpoints
            ax: a matplotlib axis
    '''
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
