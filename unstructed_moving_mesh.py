from math import gamma
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import random

#--------------------------------------------------------------------
# parameters
#--------------------------------------------------------------------
# constants
L = 1
N = 100
cfl = 0.8
gamma = 5.0/3.0
cs       = 1.0
end_time = 5.0

# plotting parameters
nstep_per_image = 1

#--------------------------------------------------------------------
# main
#--------------------------------------------------------------------
# set initial condition
t = 0
points = np.empty( (N,2) )
U = np.empty( (N,4) )
F = np.empty( (N,2,4) )

def init():
    cell = 200
    random.seed(10)
    sigma = 0.1
    listx = []
    listy = []
    while(len(listx)<N):
        for x in range(cell -1):
            for y in range(cell -1):
                if (random.random() < np.e**-(((x-cell/2)*L/cell/sigma)**2+((y-cell/2)*L/cell/sigma)**2)):
                    listx.append((x+1)/cell)
                    listy.append((y+1)/cell)

    index = np.arange(0,len(listx))
    random.shuffle(index)
    U[:,0] = 1
    U[:,1:3] = 0
    U[:,3] = 5
    for i in range(N):
        points[i,0] = listx[index[i]]
        points[i,1] = listy[index[i]]
        F[i,:] = Conserved2Flux(U[i,:])

def Conserved2Flux(u):
    flux = np.empty( (2,4) )

    P = ComputePressure( u[0], u[1], u[2], u[3] )
    wx = u[1] / u[0]
    wy = u[2] / u[0]
    flux[0,0] = u[1]
    flux[1,0] = u[2]
    flux[0,1] = wx*u[1] + P
    flux[1,2] = wy*u[2] + P
    flux[0,2] = u[1]*wy
    flux[1,1] = u[2]*wx
    flux[0,3] = (u[3]+P)*wx
    flux[1,3] = (u[3]+P)*wy

    return flux

def ComputePressure( rho, rhovx, rhovy, E ):
   P = (gamma-1.0)*( E - 0.5*(rhovx**2.0+rhovy**2)/rho )
   return P

# o------------------------------
def ComputeTimestep( U ):
   P = ComputePressure( U[:,0], U[:,1], U[:,2] )
   a = ( gamma*P/U[:,0] )**0.5
   u = np.abs( U[:,1]/U[:,0] )

   max_info_speed = np.amax( u + a )
   dt_cfl         = cfl*min_distance()/max_info_speed
   dt_end         = end_time - t

   return min( dt_cfl, dt_end )

def min_distance():
    min = L*2**0.5
    for i in range(N):
        for j in range(N):
            dis = ((points[i,0]-points[j,0])**2+ (points[i,1]-points[j,1])**2)**0.5
            if dis < min:
                min = dis
    return min
# o------------------------------

# o------------------
def ComputeNearCellLabel(points):
    nn = [ [] for i in range(N) ]
    for label in range(N):
        edgevec=vor.vertices[vor.regions[vor.point_region[label]]]
        numnode=len(edgevec)
        medge=[]
        for i in range(numnode):
            j = (i+1)%numnode
            midP = (edgevec[i]+edgevec[j])/2.
            medge.append(midP)
                
        for j in range(numnode):
            mag = np.linalg.norm(points[label]-medge[j])
            for k in range(N):
                if(k != label):
                    mag1 = np.linalg.norm(points[k]-medge[j])
                    if(np.abs(mag-mag1) < 10**-5):
                        box = mag1
                        boxlabel = k
            nn[label].append(boxlabel)
            #print(points[label],boxlabel,points[boxlabel])
    print(nn)
    return nn #如果有重複label表示不是close cell
# o------------------

# x------------------
def ComputeFluxFunction():
# x------------------

# x------------------
def ComputeW():
    return
# x------------------

# x------------------
def ComputeAveragedFluxij():
# x------------------

# x------------------
def Sum_AijFij():
# x------------------

# x------------------
def UpdatePosition():
# 0.5dt
# x------------------

# x------------------
def UpdateU():
# x------------------

# x------------------vivi
def plot(t):
    ax = plt.gca()
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    plt.savefig('fig_%03d.png'%t)
# x------------------

# o-------------------
# init point, U, F
init()
vor = Voronoi(points)
fig = voronoi_plot_2d(vor)
plot(fig, 0)
# o-------------------

num = 0
while(t<end_time):

    for i in range(nstep_per_image):
        dt = ComputeTimestep( U )
        F =  Conserved2Flux(U)
        UpdateU(dt)
        UpdatePosition(dt)
        t = dt+t
        vor = Voronoi(points)
        if(t>=end_time):break
    num +=1
    plot(num)



