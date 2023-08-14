#CAMIP geometry and rates, assuming a 1600x1600 km domain
#with the orgin at a cell center
from math import sqrt, atan2, cos, sin, pi

R_MAX = 500e+3

THULE_R = 800e3;
THULE_BC = 900.0; 
THULE_BL = -2000; 
THULE_BA = 1100;  

OFFSET=0

def translate(xb, yb):
    # translate bisicles coords (origin at cell corner)
    # to camip coords (origin at cell center)
    return xb - OFFSET, yb - OFFSET 

def thule(xb,yb, R, Bc, Bl, Ba):
    x, y = translate(xb,yb)
    rc = 0
    r = sqrt(x*x +y*y);
    theta = atan2(y, x);
    l = R - cos(2. * theta) * R/2;
    a = Bc - (Bc - Bl) * (r - rc)**2 /(R - rc)**2;
    B = Ba * cos(3*pi*r/l)+a
   
    return B
    
def bed_fullthule(x,y):
    return thule(x, y, THULE_R, THULE_BC, THULE_BL, THULE_BA)


def dome(xb,yb,R,Bc,Bl,Ba):
    x, y = translate(xb ,yb)
    rc = 0
    r = sqrt(x*x +y*y);
    B = Bc - (Bc - Bl) * (r - rc)**2 /(R - rc)**2;
    return B

def bed_fulldome(x,y):
    return dome(x, y, THULE_R, THULE_BC, THULE_BL, THULE_BA) 

def thickness_750(xb,yb):
    #constant thickness disc, radius 700 km
    x, y = translate(xb,yb)
    Xsq = x**2 + y**2
    h = 150.0 if Xsq < (750e+3)**2 else 0.0
    return h

def melt_radius_768(xb, yb, t, *etc):
    x,y = translate(xb, yb)
    Xsq = x**2 + y**2
    m = 0.0 if Xsq < (768e3)**2 else -1.0e+4
    return m

def melt_radius_750(xb, yb, t, *etc):
    x,y = translate(xb, yb)
    Xsq = x**2 + y**2
    m = 0.0 if Xsq < (750e3)**2 else -1.0e+4
    return m


def calving_rate_u_radius_750_steep(xb, yb, t, *etc):
    # calving rate cu proporional to speed
    # c = 1 on a defined radius R=750 km,
    # c << 1 for R < 750, c >> 1 for R > 750
    x, y = translate(xb, yb)
    rsq = (750e+3)**2
    z = ((x**2 + y**2)/rsq) - 1.0
    c = 0
    tol = 0.1
    c = 1.0 + z/tol
    if c < 0.0:
        c = 0.0
    elif c > 2.0: 
        c = 2.0
    return c


def calving_rate_expt34_a(xb, yb, t, *etc):
    # Want front advance w = u - c, where
    # apply rate c = au - wu/|u|
    # a = 1 for t < 1500 (fixed front for 1000 yeas, then supplied retreat)
    # a = 0 for t > 1500 (free advance) unless r > 750e+3; then 1.0
    # this function provides a
    x, y = translate(xb, yb)
    rsq = x**2 + y**2
    b = 1.0 if rsq > (750e+3)**2 else 0.0
    
    a = 1.0 if t < 1500.0 else b
    return a

def calving_rate_expt34_w(x, y, tin, thk, topg, *etc):
    # Want front advance w = u - c, where
    # w = -300 * sin (2 pi t / 1000)
    #
    # (for 1000 < t < 1500) (expt 4 retreat), and for floating ice only
    # (othersie, 0 ) (expt 3, expt 4 re-advance)
    #
    # a more physical approach:
    # apply rate c = au - wu/|u|
    #
    # this function provides -w
    # 
    
    hf = max(0, -topg * 1028.0/917.0)
    t = tin
    w = -750.0 * sin(pi * t / 500.0) if (thk < hf and t > 1000.0 and t < 1500.0) else 0.0
    return -w
