import math


ylen = 16.0e+3
xlen = 128.0e+3

def thickness(x,y,*etc):
    
    h = 0.0 if x > 96.0e+3 else 100.0
    return h 

def topography(x,y,*etc):

    x0 = 0.0
    x1 = 96.0e+3
    b = 100.0 * (x - x1)/(x0 - x1) 
    
    return b
            
def acab(x,y,t,thck,topg,*etc):
    return 0.5


t_retreat = 8.0

def marine_retreat_simpleP(x,y,t,thck,topg,*etc):

    p = 0.0
    if (topg < 0.0) and (t > t_retreat):
        p = 1.0
        
    return p

def marine_retreat_simpleI(x,y,t,thck,topg,*etc):

    retreat_rate = 0.0
    if (topg < 0.0) and (t > t_retreat):
        retreat_rate = 128.0
        
    return retreat_rate
