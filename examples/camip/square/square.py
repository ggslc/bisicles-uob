from math import sqrt

def velocity(x,y,*etc):
    vx = -500.0
    vy = 0.0
    return (vx,vy)

def thickness(x,y):
    H = 1.0e-4
    if x > 32.0e3 and y < 32.0e3:
        H = 1.0e+3
    return H

def topography(x,y):
    return 1.0

def constfriction(x,y,t,thck,topg):
    return 1.0e+4

def acab(x,y,*etc):
    return 0.0

def rate1000(x,y,t,thck,topg):
    return 1.0e+3

def rate750(x,y,t,thck,topg):
    return 750.0

def rate500(x,y,t,thck,topg):
    return 500.0

def rate250(x,y,t,thck,topg):
    return 250.0

def rate0(x,y,t,thck,topg):
    return 0.0


