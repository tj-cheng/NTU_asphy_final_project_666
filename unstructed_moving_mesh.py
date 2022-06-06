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
end_time = 1.0

# plotting parameters
nstep_per_image = 10

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
    U[:,3] = 50
    for i in range(N):
        points[i,0] = listx[index[i]]
        points[i,1] = listy[index[i]]
        F[i] = Conserved2Flux(U[i])

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
   if P<0:
    #    print("P<0")
       P = 0
   return P

# o------------------------------
def ComputeTimestep( U ):
    P = np.empty(N)
    for i in range(N):
        P[i] = ComputePressure( U[i,0], U[i,1], U[i,2], U[i,3] )
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
            if(i!=j):
                dis = ((points[i,0]-points[j,0])**2+ (points[i,1]-points[j,1])**2)**0.5
                if dis < min:
                    min = dis
    return min
# o------------------------------

# o------------------
def ComputeNearCellLabel(nth_point):
    fp = [] # face point
    edgevec=vor.regions[vor.point_region[nth_point]]
    numnode=len(edgevec)
    medge = []
    for i in range(numnode):
        j = (i+1)%numnode
        if(edgevec[i]!=-1 and edgevec[j]!=-1):
            midP = (vor.vertices[edgevec[i]]+vor.vertices[edgevec[j]])/2.
            mag = np.linalg.norm(points[nth_point]-midP) # point to midpoint distance
            medge.append(midP)
            for k in range(N):
                if(k != nth_point):
                    mag1 = np.linalg.norm(points[k]-midP)
                    if(np.abs(mag-mag1) < 1e-15):
                        fp.append(k)
    if len(fp)!= len(medge): print("error", nth_point, len(fp), len(medge))
    return fp, medge
# o------------------

# o------------------
def ComputeFluxFunction(nth_point, face_point,dt):
    dx = points[nth_point]-points[face_point]
    fflux = 0.5*(F[nth_point]+F[face_point]-(U[nth_point]-U[face_point])*np.transpose(np.array([dx]))/dt)
    return fflux
# o------------------

# o------------------
def ComputeW(nth_point, face_point, MidEdge):
    w_this = (U[nth_point,1:3] / U[nth_point,0])
    w_face = (U[face_point,1:3] / U[face_point,0])
    r_this = points[nth_point]
    r_face = points[face_point]
    
    wp = (w_this-w_face)*(MidEdge-(r_this+r_face)/2)*(r_face-r_this)/(np.linalg.norm(r_this-r_face)**2)
    w = (w_this+w_face)/2 +wp
    return w
# o------------------

# x------------------
def ComputeAveragedFlux(nth_point, mth_edge, face_point, MidEdge, dt):
    f = ComputeFluxFunction(nth_point, face_point, dt)
    w = ComputeW(nth_point, face_point, MidEdge)
    Aij_vec = ((vor.vertices[vor.regions[vor.point_region[nth_point]][mth_edge]]-MidEdge)*2).dot(np.array([[0,1],[-1,0]]))
    Fij = np.array([Aij_vec]).dot((f-U[nth_point]*np.transpose(np.array([w]))))
    Fij = Fij.flatten()
    return Fij

# x------------------

# o------------------
def UpdatePosition(dt):
    for i in range(N):
        wx = (U[i,1] / U[i,0])
        wy = (U[i,2] / U[i,0])
        points[i,0] += wx*dt
        points[i,1] += wy*dt
# o------------------

# o------------------
def UpdateU(dt):
    for i in range(N):
        face_point, MidEdge = ComputeNearCellLabel(i)
        # print(len(face_point))
        for j in range(len(face_point)):
            flux = ComputeAveragedFlux(i, j, face_point[j], MidEdge[j], dt)
            U[i] -= flux*dt
# o------------------

# x------------------vivi
def plot(vor,t):
    fig = voronoi_plot_2d(vor,show_vertices = False)
    ax = plt.gca()
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_aspect('equal', adjustable='box')
    plt.savefig('fig_%03d.png'%t,dpi=150)
# x------------------

# o-------------------
# init point, U, F
init()
vor = Voronoi(points)
plot(vor, 0)
# o-------------------
dt = ComputeTimestep( U )
num = 0
while(t<end_time):

    for i in range(nstep_per_image):
        dt = ComputeTimestep( U )
        for j in range(N):
            F[j] =  Conserved2Flux(U[j])
        UpdateU(0.5*dt)
        UpdatePosition(dt)
        vor = Voronoi(points)
        for j in range(N):
            F[j] =  Conserved2Flux(U[j])
        UpdateU(0.5*dt)
        t = dt+t
        vor = Voronoi(points)
        print("%e -> %e , dt = %e"%(t-dt, t, dt))
        if(t>=end_time):break
    num +=1
    plot(vor, num)

# for i in range(1):

#     dt = ComputeTimestep( U )
#     for j in range(N):
#         F[j] =  Conserved2Flux(U[j])
#     UpdateU(dt)
#     print(U[0])
#     UpdatePosition(dt)
#     t = dt+t
#     vor = Voronoi(points)
#     print("%e -> %e , dt = %e"%(t-dt, t, dt))



