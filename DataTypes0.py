#-------------------------------- DataTypes0.py -------------------------------
#  This file is part of the MDsquad simulation package.
#  Defines the main Atom class
#  Copyright Titus Beu, 2019
#------------------------------------------------------------------------------
from copy import *
from Vectors import *

class Atom:
   def __init__(self, r = VecR3(),
                name = "", type = "", resi = "", chrg = 0e0):
      self.name = name                                              # atom name
      self.symb = name[0] if len(name) else ""                # chemical symbol
      self.type = type                                              # atom type
      self.resi = resi                                           # residue name
      self.ires = 1                                      # residue sequence no.
      self.segm = ""                                             # segment name
      self.mass = 0e0                                         # atom mass (amu)
      self.chrg = chrg                                        # atom charge (e)
      self.occp = 0e0                   #  occupancy - user field usable in VMD
      self.beta = 0e0                                    #  temeperature factor
      self.r = deepcopy(r)                 # position (A) (field-by-field copy)

class Bond:
   def __init__(self, indx = (0,0)):
      self.indx = indx                                           # atom indexes
