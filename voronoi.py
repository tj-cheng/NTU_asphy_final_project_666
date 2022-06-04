import numpy as np
import matplotlib
import matplotlib.pyplot as plt

points = np.array([[3, .2], [-2, 1], [0, 2.1], [1.1, -1], [1, 1], [0.9, 1.9],[2.2, -0.2], [1.9, 1.1], [2.1, 2.3]])

from scipy.spatial import Voronoi, voronoi_plot_2d
vor = Voronoi(points)
label = 4

print("points[4] =", points[label])
#print("ridge_points=",vor.ridge_points[5])
#print("ridge_vertics=",vor.ridge_vertices[5])
print("region for point 4 =", vor.point_region[label])
print("vertices for region for point 4 =", vor.regions[vor.point_region[label]])
print(points[vor.regions[vor.point_region[label]]])
print("coordinate of Voronoi vertices associated with point 4:")
print(vor.vertices[vor.regions[vor.point_region[label]]])
edgevec=vor.vertices[vor.regions[vor.point_region[label]]]
numnode=len(vor.vertices[vor.regions[vor.point_region[label]]])
mag=0.
mag1=0.
box = 10000
boxlabel=0
medge=[]
for i in range(numnode):
    j = (i+1)%numnode
    midP = (edgevec[i]+edgevec[j])/2.
    medge.append(midP)
print(medge)

        
for j in range(numnode):
    mag = np.linalg.norm(points[label]-medge[j])
    box = 10000
    for k in range(9):
        if(k != label):
            mag1 = np.linalg.norm(points[k]-medge[j])
            if(np.abs(mag-mag1) < 10**-5):
                box = mag1
                boxlabel = k
    print(boxlabel,points[boxlabel])



        
        
import matplotlib.pyplot as plt
fig = voronoi_plot_2d(vor)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.show()

