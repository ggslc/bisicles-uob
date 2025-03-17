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


t_retreat = 8.0

def marine_retreat_simpleP(x,y,t,thck,topg,*etc):
    
    return 1.0

def marine_retreat_simpleI(x,y,t,thck,topg,*etc):

    retreat_rate = 0.0
    if (topg < 0.0) and (t > t_retreat):
        retreat_rate = 128.0
        
    return retreat_rate

def marine_retreat_simpleN(x,y,t,thck,topg,*etc):

    retreat_rate = 0.0
    if (topg < 0.0) and (t > t_retreat):
        retreat_rate = 128.0
        
    return retreat_rate
