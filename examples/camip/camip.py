#CAMIP geometry and rates
from math import sqrt, atan2, cos, sin, pi

R_MAX = 500e+3

THULE_R = 800e3;
THULE_BC = 900.0; 
THULE_BL = -2000; 
THULE_BA = 1100;  

JIM_R = 400e3;
JIM_BC = 450.0; 
JIM_BL = -1000; 
JIM_BA = 550; 

def thule(x,y, R, Bc, Bl, Ba):

    rc = 0
    r = sqrt(x*x +y*y);
    theta = atan2(y, x);
    l = R - cos(2. * theta) * R/2;
    a = Bc - (Bc - Bl) * (r - rc)**2 /(R - rc)**2;
    B = Ba * cos(3*pi*r/l)+a
   
    return B
    
def bed_thule_hilmar(x,y):
    return thule(x, y, THULE_R, THULE_BC, THULE_BL, THULE_BA)

def bed_thule_jim(x,y):
    return thule(x, y, JIM_R, JIM_BC, JIM_BL, JIM_BA)

def dome(x,y,R,Bc,Bl,Ba):
    rc = 0
    r = sqrt(x*x +y*y);
    theta = atan2(y, x);
    B = Bc - (Bc - Bl) * (r - rc)**2 /(R - rc)**2;
    return B

def bed_dome(x,y):
    return dome(x, y, THULE_R, THULE_BC, THULE_BL, THULE_BA) 

def thickness_768(x,y):
    #constant thickness disc, radius 768 km
    Xsq = x**2 + y**2
    h = 250.0 if Xsq < (768e+3)**2 else 0.0
    return h



def thickness_500(x,y):
    #constant thickness disc, radius 500 km
    Xsq = x**2 + y**2
    h = 250.0 if Xsq < (500e+3)**2 else 0.0
    return h

def thickness_600(x,y):
    #constant thickness disc, radius 600 km
    Xsq = x**2 + y**2
    h = 150.0 if Xsq < (600e+3)**2 else 0.0
    return h

def topography_dome_0(x,y):
    #my dome - need Jim's
    Xsq = x**2 + y**2
    Rsq = (1000e+3)**2
    z = -300 + 400.0 * (1.0 - Xsq/Rsq)**0.5
    return z


def calving_rate_0(x,y,t,*etc):
    cr = 1000.0
    return 0.0

def melt_radius_750(x, y, t, *etc):
    Xsq = x**2 + y**2
    m = 0.0 if Xsq < (750e3)**2 else -1.0e+4
    return m

def mild_radius_750(x, y, t, *etc):
    Xsq = x**2 + y**2
    m = 0.0 if Xsq < (750e3)**2 else -0.3e4
    return m


def melt_radius_768(x, y, t, *etc):
    Xsq = x**2 + y**2
    m = 0.0 if Xsq < (768e3)**2 else -1.0e+4
    return m

def calving_rate_u1p1(x, y, t, *etc):
    # calving rate cu = 1.1u
    return 3.0

def calving_rate_zero(x, y, t, *etc):
    # calving rate cu = 0u
    return 0.0


def calving_rate_u_radius_750(x, y, t, *etc):
    # calving rate cu proporional to speed
    # c = 1 on a defined radius R=750 km
    # > 1 for R > 750, < 1 for R < 750 lm
    rsq = (750e+3)**2
    c = ((x**2 + y**2)/rsq)
    return c

def calving_rate_u_radius_750_steep(x, y, t, *etc):
    # calving rate cu proporional to speed
    # c = 1 on a defined radius R=750 km,
    # c << 1 for R < 750, c >> 1 for R > 750
    rsq = (750e+3)**2
    z = ((x**2 + y**2)/rsq) - 1.0
    c = 0
    tol = 0.25
    c = 1.0 + z/tol
    if c < 0.0:
        c = 0.0
    elif c > 2.0: 
        c = 2.0
    return c

def calving_rate_expt2_brutal(x, y, t, *etc):
    # Want front advance w = u - c, where
    # w = -300 * sin (2 pi t / 1000)
    #
    # 
    # This is a brute-force approach: move
    # the expected front position
    # r = r0 + (A cos (omega t) - 1)
    # (which satifies dr/dt = w and r(0) = r0)
    # then apply a calving rate cu > u at points
    # outside the circle with radius r, 
    # and cu < u inside (and cu = u right on it)
    
    r0 = 750e+3
    omega = 2.0 * pi / 1.0e+3
    amp = 300.0 / omega
    r = r0 + amp * (cos (omega * t) - 1.0)
    z = sqrt(x**2 + y**2)/r - 1.0
    c = 0
    tol = 1.0e-2
    c = 1.0 + z/tol
    if c < 0.0:
        c = 0.0
    elif c > 2.0: 
        c = 2.0
    return c

def calving_rate_expt2_a(x, y, t, *etc):
    # Want front advance w = u - c, where
    # w = 300 * sin (2 pi t / 1000)
    #
    # a more physical approach:
    # apply rate c = au - wu/|u|
    #
    # this function provides a    
    return 1.0

def calving_rate_expt2_w(x, y, t, *etc):
    # Want front advance w = u - c, where
    # w = -300 * sin (2 pi t / 1000)
    #
    # a more physical approach:
    # apply rate c = au - wu/|u|
    #
    # this function provides -w
    w = -300.0 * sin(pi * t / 500.0) 
    return -w
