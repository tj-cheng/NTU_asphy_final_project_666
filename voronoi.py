import numpy as np
#points = np.array([[0, .2], [-.2, 1], [0, 2.1], [1.1, -.1], [1, 1], [0.9, 1.9],[2.2, -0.2], [1.9, 1.1], [2.1, 2.3]])
points = np.array([[3, .2], [-.2, 1], [0, 2.1], [1.1, -.1], [1, 1], [0.9, 1.9],[2.2, -0.2], [1.9, 1.1], [2.1, 2.3]])

from scipy.spatial import Voronoi, voronoi_plot_2d
vor = Voronoi(points)

print("points[5] =", points[5])
print(vor.nlist[points[5]])
print("region for point 5 =", vor.point_region[5])
print("vertices for region for point 5 =", vor.regions[vor.point_region[5]])
print("coordinate of Voronoi vertices associated with point 5:")
print(vor.vertices[vor.regions[vor.point_region[5]]])

import matplotlib.pyplot as plt
fig = voronoi_plot_2d(vor)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.show()
