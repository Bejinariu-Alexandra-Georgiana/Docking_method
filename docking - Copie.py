from Tool import *
from grafxlib import *

GraphInit(800,800)
a = VecR3(3, 4, 8)
b = VecR3(6, 7, 9)
b1 = VecR3(6, 7, 9)
c = VecR3(8.6, 4.5, 6.5)
c1 = VecR3(8.6, 4.5, 6.5)

(d, e, f, p1, p2, p0, C1) = Tools0(a, b, c)


j = VecR3(1, 2, 3)
k = VecR3(6 , 4, 1)
k1 = VecR3(6 , 4, 1)
l = VecR3(3 , 5.196, 0)
l1 = VecR3(3 , 5.196, 0)

(g, h, i, p3, p4, p5, C2) = Tools0(j, k, l)


d = 10
g.z = g.z + d
h.z = h.z + d
i.z = i.z + d

n = 6
x = [0, a.x, b1.x, c1.x, j.x, k1.x, l1.x]
y = [0, a.y, b1.y, c1.y, j.y, k1.y, l1.y]
z = [0, a.z, b1.z, c1.z, j.z, k1.z, l1.z]

n3 = 2
ind1 = [0, 1, 4]
ind2 = [0, 2, 5]
ind3 = [0, 3, 6]


title = "initial "
PlotStruct(x,y,z,n,ind1,ind2,ind3,n3,0e0,0e0,0,0e0,0e0,0,0.05,0.45,0.55,0.92,title)

n = 6
x = [0, a.x, b.x, c.x, j.x, k.x, l.x]
y = [0, a.y, b.y, c.y, j.y, k.y, l.y]
z = [0, a.z, b.z, c.z, j.z, k.z, l.z]

n3 = 2
ind1 = [0, 1, 4]
ind2 = [0, 2, 5]
ind3 = [0, 3, 6]

title = "final "

PlotStruct(x,y,z,n,ind1,ind2,ind3,n3,
           0e0,0e0,0,0e0,0e0,0,0.55,0.95,0.05,0.42,title)
