#--------------------------------- Vectors.py ---------------------------------
#  This file is part of the MDsquad simulation package.
#  Defines an R3 vector structure and applicable vector operations.
#  Copyright Titus Beu, 2017
#------------------------------------------------------------------------------
from math import *

class VecR3:
   def __init__(self, x = 0e0, y = 0e0, z = 0e0):
      self.x = x; self.y = y; self.z = z                 # Cartesian components

def VecZero():                                                    # Zero vector
   return VecR3(0e0, 0e0, 0e0)

def VecFac(f, v):                                   # Multiply vector by scalar
   w = VecR3()
   w.x = f * v.x
   w.y = f * v.y
   w.z = f * v.z
   return w

def VecDiv(v, f):                                     # Divide vector by scalar
   w = VecR3()
   w.x = v.x / f
   w.y = v.y / f
   w.z = v.z / f
   return w

def VecSum(u, v):                                    # Sum w of vectors u and v
   w = VecR3()
   w.x = u.x + v.x
   w.y = u.y + v.y
   w.z = u.z + v.z
   return w

def VecDif(u, v):                             # Difference w of vectors u and v
   w = VecR3()
   w.x = u.x - v.x
   w.y = u.y - v.y
   w.z = u.z - v.z
   return w

def VecLin1(u, f, v):                            # Linear combination u + f * v
   w = VecR3()
   w.x = u.x + f * v.x
   w.y = u.y + f * v.y
   w.z = u.z + f * v.z
   return w

def VecLin2(f, u, g, v):                     # Linear combination f * u + g * v
   w = VecR3()
   w.x = f * u.x + g * v.x
   w.y = f * u.y + g * v.y
   w.z = f * u.z + g * v.z
   return w

def VecDot(u, v):                                                 # Dot product
   return (u.x * v.x + u.y * v.y + u.z * v.z)

def VecNorm(v):                                                # Norm of vector
   return sqrt(v.x * v.x + v.y * v.y + v.z * v.z)

def VecUnit(v):                                                   # Unit vector
   n = VecNorm(v)
   return VecFac((1e0/n if n else 0e0), v)

def VecCos(u, v):                             # Cosine of angle between vectors
   n = VecNorm(u) * VecNorm(v)
   return (VecDot(u,v) / n if n else 0e0)

def VecCross(u, v):                                             # Cross product
   w = VecR3()
   w.x = u.y * v.z - u.z * v.y
   w.y = u.z * v.x - u.x * v.z
   w.z = u.x * v.y - u.y * v.x
   return w

def VecSin(u, v):                               # Sine of angle between vectors
   n = VecNorm(u) * VecNorm(v)
   return (VecNorm(VecCross(u,v)) / n if n else 0e0)

# Wrap vector r within [-Lbox/2,+Lbox/2]
def VecWrap(r, Lbox, Lbox2):
   if   (r.x >  Lbox2.x): r.x -= Lbox.x
   elif (r.x < -Lbox2.x): r.x += Lbox.x
   if   (r.y >  Lbox2.y): r.y -= Lbox.y
   elif (r.y < -Lbox2.y): r.y += Lbox.y
   if   (r.z >  Lbox2.z): r.z -= Lbox.z
   elif (r.z < -Lbox2.z): r.z += Lbox.z
