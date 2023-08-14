from math import sqrt

L = 640.0e+3 #domain half width/length

rf = 256.0e+3
H = 1000.0
D = -500.0
umax = 500.0

def get_r(x,y):
    return sqrt((x-L)**2+(y-L)**2)

def velocity(x,y,*etc):
    r_cos_theta = x - L
    r_sin_theta = y - L
    vx = umax/rf * r_cos_theta
    vy = umax/rf * r_sin_theta
    return (vx,vy)

def thickness(x,y):
    if (x-L)**2 + (y-L)**2 > rf**2:
        return 0.0
    return H

def topography(x,y):
    return D

def constfriction(x,y,t,thck,topg):
    return 1.0e+4

def acab(x,y,*etc):
    r = get_r(x,y)
    a = 2*H*umax/rf
    #if ((x-L)**2 + (y-L)**2) > rf**2:
    #    return -1000000.0
    return a

def rate1000(x,y,t,thck,topg):
    return 1.0e+3

def rate750(x,y,t,thck,topg):
    return 750.0

def irate750(x,y,t,thck,topg):
    u,v = velocity(x,y,t)
    uu = 1.5 * (u*u + v*v)**0.5
    return uu


def rate500(x,y,t,thck,topg):
    return 500.0

def rate250(x,y,t,thck,topg):
    return 250.0

def rate0(x,y,t,thck,topg):
    return 0.0


