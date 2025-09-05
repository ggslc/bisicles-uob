import math as m

L = 128.0e+3
W = 128.0e+3
H0 = 512.0 # thickness at left boundary
HL = 256.0 # thickness at CF
B0 = -HL # bed elevation at left boundary
BL = -H0 # bed elevation at CF

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
    h = H0*(x-L)/(-L) + HL*x/L
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


