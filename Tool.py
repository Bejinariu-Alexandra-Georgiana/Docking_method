from Vectors import *
from Tools_lib import *
from DataTypes0 import *
from DataInOut0 import *

def Tools0(a = VecR3(0e0, 0e0, 0e0), b = VecR3(0e0, 0e0, 0e0), c = VecR3(0e0, 0e0, 0e0)):

    ox = VecR3(1e0,0e0,0e0)                          # 3 axis-->ox, oy, oz
    oy = VecR3(0e0,1e0,0e0)
    oz = VecR3(0e0,0e0,1e0)

    w = VecDif(b,a)                                               # vectors 
    v = VecDif(c,a)


    n = VecCross(w,v)                                         # normal vector of the plane


    q = VecCross(n, oz)                                         # normal vector between normal vector of the plane and axis z 


    u1 = UnitV(q)                                           # unit vector of q

    teta = acos(VecCos(n, oz))                                     # the angle between z si normal vector of the plane

    matrix1 = Mrot(teta, u1)                                   

    nf = Rotation(n, matrix1)                                   # nf || z

    vf = Rotation(v, matrix1)                                   # vf || xoy                                  


    wf = Rotation(w, matrix1)                                 # wf ||xoy


    r = VecSum(vf, wf)

    alfa = acos(VecCos(r, ox))
    matrix2 = Mrot(alfa, UnitV(VecCross(r, ox)))                 
    rf = Rotation(r, matrix2)                                # rf || ox


    vf1 = Rotation(vf, matrix2)

    wf1 = Rotation(wf, matrix2)

    a = VecZero()

    b.x = vf1.x
    b.y = vf1.y
    b.z = vf1.z

    c.x = wf1.x
    c.y = wf1.y 
    c.z = wf1.z

    C = findCircle(a, b, c)

    #print("C initial=", C.x, C.y, C.z)
    
    r = VecZero()
    MoveC(a, C, r)
    MoveC(b, C, r)
    MoveC(c, C, r)

    C = findCircle(a, b, c)
    #print("C final = ", C.x, C.y, C.z)

    wf1 = VecDif(b, a)
    vf1 = VecDif(c, a)
    rf1 = VecSum(vf1, wf1)
    
    return a, b, c, vf1, wf1, rf1, C















