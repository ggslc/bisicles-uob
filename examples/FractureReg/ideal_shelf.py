import math as m

L = 128.0e+3
W = 128.0e+3
S0 = 512.0 # elevation at left boundary
SL = 256.0 # elevation at CF
B0 = -SL # bed elevation at left boundary
BL = -S0 # bed elevation at CF

# crack parallel to front
L3 = 32.0e+3 #crack to front
L2 = 4.0e+3  #crack width
L1 = L - L3 - L2  #left hand boundary to crack

#crack parallel to shear margin
W1 = 32.0e+3 # y = 0 to crack
W2 = 4.0e+3  # crack width
W3 = W - W1 + W2


def topography(x,y):
    return B0*(x - L)/(-L) + BL*x/L

def thickness_basic(x,y):
    #basic wedge shape
    s = S0*(x-L)/(-L) + SL*x/L
    b = topography(x,y)
    h = min(s-b,s)
    #strip if ocean
    if (x > L - 2.0e+3):
        h = 0.0
    return h
   
   
def mucoef_rifted(x,y,*etc):

    mucoef = 1.0

    #rift parallel to front
    if (x > L1) and (x < L1 + L2) and (y > W1) and (y < W - W1):
        mucoef = 0.25
            
    #(higher) damage parallel to shear margin @ y = W1
    if (y > W1) and (y < W1 + W2):
        mucoef = min(mucoef,0.25)
    #(lower) damage parallel to shear margin @y = W - W1
    if (y > W - W1 - W2) and (y < W - W1):
        mucoef = min(mucoef, 0.5)

    return mucoef

BETA_MAX = 2.0e+3
BETA_MID = 1.0e+3

def beta(x,y,*etc):

    beta = BETA_MAX
    if (y > W1) and (y < W - W1):
        beta = BETA_MID * (1.0 + m.cos(16.0 * m.pi * x/L))

    return beta

DIRT = 10.0 # amplitude of faulty 
H_EPS = 50.0 # low thickness

def preprocess_synthetic_data(s, b, h, ux, uy, x, y, *etc):

    #provide inverse problem data given 'observations'
    #surface elevation (s)
    #bedrock elevation (b)
    #thickness (h)
    #velocity (ux, uy)
    #co-ordinates (x,y)

    beta_init = BETA_MID 
    umod_clean = m.sqrt((ux*ux + uy*uy))
    yy = m.pi * (x + y) * 32.0/W
    umod_dirty = umod_clean + (DIRT * (m.sin(yy)) if s - h > b else 0.0)
    uc_grounded = 0.0 if s - h > b + H_EPS else 1.0
    uc_all = 0.0 if h < H_EPS else 1.0

    return beta_init, umod_clean, umod_dirty, uc_grounded, uc_all


