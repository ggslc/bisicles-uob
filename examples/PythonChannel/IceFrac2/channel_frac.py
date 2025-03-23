import math

ela = 100.0
ylen = 16.0e+3
xlen = 128.0e+3
rhoi = 917.0
rhow = 1028.0
frat = (1-rhoi/rhow)

def thickness(x,y,*etc):
    
    h = 0.0 if x > 96.0e+3 else 100.0
    return h 

def topography(x,y,*etc):

    x0 = 0.0
    x1 = 72.0e+3
    b = 100.0 * (x - x1)/(x0 - x1) 
    
    return b
            
def acab(x,y,t,thck,topg,*etc):
    s = max(thck+topg, frat*thck)
    return 0.5/100.0 * (s - ela)


t_advance = 8.0
t_hold = 16.0
t_retreat = 24.0 # a
r_retreat = 128.0 # m/a
u_typical = 256.0 # m/a

def marine_retreat(topg, t):
    retreat_rate = 0.0
    if (topg < 0.0) and (t > t_retreat):
        retreat_rate = r_retreat
    return retreat_rate
    
def hold_advance_hold(t):
    uc = 1.0
    if t > t_advance:
        uc = 0.0
    if t > t_hold:
        uc = 1.0
        
    return uc

def hold(t):
    return 1.0



#'rn' experiments - retreat normal to front
def marine_retreat_simplePrn(x,y,t,thck,topg,*etc):
    return hold(t)
def marine_retreat_simpleIrn(x,y,t,thck,topg,*etc):        
    return 0.0
def marine_retreat_simpleNrn(x,y,t,thck,topg,*etc):
    return marine_retreat(topg, t)

#'ri' experiments - retreat opposed to velocity, set magnitude
def marine_retreat_simplePri(x,y,t,thck,topg,*etc):
    return hold(t)
def marine_retreat_simpleIri(x,y,t,thck,topg,*etc):        
    return marine_retreat(topg,t)
def marine_retreat_simpleNri(x,y,t,thck,topg,*etc):
    return 0.0

#'rp' experiments - retreat opposed to velocity, proporional magnitude
def marine_retreat_simplePrp(x,y,t,thck,topg,*etc):
    prop = hold(t) + marine_retreat(topg,t)/u_typical    
    return prop
def marine_retreat_simpleIrp(x,y,t,thck,topg,*etc):        
    return 0.0
def marine_retreat_simpleNrp(x,y,t,thck,topg,*etc):
    return 0.0



