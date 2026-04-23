C3_MAX = 6.0e4
C3_SEA = 1.0e+2

TB_MAX = 1.0e+5 # 100 Kpa
UU_MAX = 1.0e+4 # stupid speed...

C2_MAX = 4.0e+4
C4_MAX = 1.0e+5
C2_SEA = C3_SEA
C4_SEA = C3_SEA

VS = [(2,C2_MAX,C2_SEA),(3,C3_MAX,C3_SEA),(4,C4_MAX,C4_SEA)]

def cm(tb, uu, h, n, c_max, c_sea):
    cm = tb/(uu**(1/n) + 1.0e-1) if h > 0.0 else c_sea
    cm = cm if uu < UU_MAX else c_max
    cm = min(cm,c_max)
    return cm

def pp(c1,ux,uy,h,*etc):

    uu = (ux**2 + uy**2)**0.5
    tb = uu*c1 if uu > 2.0 else TB_MAX
    tb = min(tb, TB_MAX)

    c2,c3,c4 = [cm(tb, uu, h, *v) for v in VS]

    return c1,c2,c3,c4,uu,tb
