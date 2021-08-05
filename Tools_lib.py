from Vectors import *
from math import *

def findCircle( a = VecZero(), b = VecZero(), c = VecZero()):
   xab = a.x - b.x
   xac = a.x - c.x

   yab = a.y - b.y
   yac = a.y - c.y

   yca = c.y - a.y
   yba = b.y - a.y

   xca = c.x - a.x
   xba = b.x - a.x

   sxac = a.x*a.x - c.x*c.x
   syac = a.y*a.y - c.y*c.y

   sxba = b.x*b.x - a.x*a.x
   syba = b.y*b.y - a.y*a.y

   f = ( sxac * xab + syac * xab + sxba * xac + syba * xac)/( 2 *(yca * xab - yba * xac))
   g = ( sxac * yab + syac * yab + sxba * yac + syba * yac)/ ( 2 *( xca * yab - xba *yac))

   C = VecZero()
   C.x = -g
   C.y = -f
   C.z = 0.0

   return C


   
def MoveC(self, v = VecR3(0e0,0e0,0e0), r = VecR3(0e0,0e0,0e0)):

   self.x += r.x - v.x
   self.y += r.y - v.y
   self.z += r.z - v.z

def UnitV(self):

   norma = VecNorm(self)

   if(norma != 0):
      u = VecDiv(self, norma)
      return u
   else:
      u = VecR3(0e0, -1e0, 0e0)

      return u
      
def Mrot(angle, u = VecR3(0e0,0e0,0e0)):
   Rot = [[0]*3 for i in range(3)]

   Rot[0][0] = cos(angle) + u.x*u.x*(1-cos(angle))
   Rot[0][1] = u.x*u.y*(1-cos(angle)) - u.z*sin(angle)
   Rot[0][2] = u.x*u.z*(1-cos(angle)) + u.y*sin(angle)
   Rot[1][0] = u.x*u.y*(1-cos(angle)) + u.z*sin(angle)
   Rot[1][1] = cos(angle) + u.y*u.y*(1-cos(angle))
   Rot[1][2] = u.y*u.z*(1-cos(angle)) - u.x*sin(angle)
   Rot[2][0] = u.z*u.x*(1-cos(angle)) - u.y*sin(angle)
   Rot[2][1] = u.z*u.y*(1-cos(angle)) + u.x*sin(angle)
   Rot[2][2] = cos(angle) + u.z*u.z*(1-cos(angle))


   return Rot

def Rotation(self, matrix):
   nf = [[0] for i in  range (3)]       #normala rotata

   ni = [[0] for i in range (3)]        #normala initiala

   ni[0][0] = self.x
   ni[1][0] = self.y
   ni[2][0] = self.z

   for i in range (3):
      for j in range(1):
         for k in range(3):
            nf[i][j] += matrix[i][k] * ni[k][j]

   nm = VecR3(0e0,0e0,0e0)
   nm.x = nf[0][0]
   nm.y = nf[1][0]
   nm.z = nf[2][0]

   return nm


   
   
   
   
   
