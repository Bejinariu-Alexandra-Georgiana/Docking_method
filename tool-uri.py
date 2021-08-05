from Vectors import *
from Tools_lib import *
from DataTypes0 import *
from DataInOut0 import *


a = VecR3(0, 3.25, 5)                      # 3 points
b = VecR3(0, 0.5, 4)
c = VecR3(0, 10.5, 4)

ox = VecR3(1e0,0e0,0e0)                          # 3 axis-->ox, oy, oz
oy = VecR3(0e0,1e0,0e0)
oz = VecR3(0e0,0e0,1e0)

file = "Triangles.pdb"
atoms = [Atom()]

atoms.append(Atom(VecR3(a.x,a.y,a.z), "C"))
atoms.append(Atom(VecR3(b.x,b.y,b.z), "C"))
atoms.append(Atom(VecR3(c.x,c.y,c.z), "C"))

w = VecDif(b,a)                                                # vectors 
v = VecDif(c,a)


n = VecCross(w,v)                                         # normal vector of the plane
print("normala initiala=", n.x, n.y, n.z)

q = VecCross(n,oz)                                       # normal vector between normal vector of the plane and axis z 
print("normala dintre normala si z=", q.x, q.y, q.z)

u1 = UnitV(q)                                           # unit vector of q
print("unit vector", u1.x, u1.y, u1.z)

teta = acos(VecCos(n, oz))                                     # the angle between z si normal vector of the plane
print("teta=", teta)

matrix1 = Mrot(teta, u1)                                   
print(matrix1)

nf = Rotation(n, matrix1)                                   # nf || z
print("nf=", nf.x, nf.y, nf.z)

vf = Rotation(v, matrix1)                                   # vf || xoy                                  
print("vf =",vf.x, vf.y, vf.z)


print(w.x , w.y, w.z)      
wf = Rotation(w, matrix1)                                 # wf ||xoy
print("wf =",wf.x, wf.y, wf.z)


r = VecSum(vf, wf)
print("rezultanta= ", r.x, r.y, r.z)

alfa = acos(VecCos(r, ox))
matrix2 = Mrot(alfa, UnitV(VecCross(r, ox)))                 
rf = Rotation(r, matrix2)                                # rf || ox

print("rf =", rf.x, rf.y, rf.z)

vf1 = Rotation(vf, matrix2)
print("vf1=", vf1.x, vf1.y, vf1.z)

wf1 = Rotation(wf, matrix2)
print("wf1=", wf1.x, wf1.y, wf1.z)

a = VecZero()

b.x = vf1.x
b.y = vf1.y
b.z = vf1.z

c.x = wf1.x
c.y = wf1.y
c.z = wf1.z

C = findCircle(a, b, c)


print("C initial=",C.x, C.y, C.z)

r = VecZero()         
MoveC(a, C, r)                                  # move the center in the origin 
MoveC(b, C, r)
MoveC(c, C, r)

C = findCircle(a, b, c)
print("C final=", C.x, C.y, C.z)

wf1 = VecDif(b, a) 
vf1 = VecDif(c, a)

print("a=", a.x, a.y, a.z)
print("b =", b.x, b.y , b.z)
print("c= ", c.x , c.y , c.z)


print("finalul wf1=", wf1.x, wf1.y, wf1.z)
print("finalul vf1=", vf1.x, vf1.y, vf1.z)

rf1 = VecSum(vf1, wf1)
print("rezultanta finala=", rf1.x, rf1.y, rf.z)


atoms.append(Atom(VecR3(a.x,a.y,a.z), "S"))
atoms.append(Atom(VecR3(b.x,b.y,b.z), "S"))
atoms.append(Atom(VecR3(c.x,c.y,c.z), "S"))




WritePDB(file, atoms)


















