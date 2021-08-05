#-------------------------------- DataInOut0.py -------------------------------
#  This file is part of the MDsquad simulation package.
#  Contains routines for reading/writing from/to pdb files
#  Copyright Titus Beu, 2019
#------------------------------------------------------------------------------
from struct import *    
from DataTypes0 import *

#==============================================================================
def ReadPDB(file):
#------------------------------------------------------------------------------
#  Reads a pdb file
#  http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
#   1 -  6     Record name   "ATOM  "
#   7 - 11     Integer       serial       Atom serial number
#  13 - 16     Atom          name         Atom name
#  17          Character     altLoc       Alternate location indicator
#  18 - 20     Residue name  resName      Residue name
#  22          Character     chainID      Chain identifier
#  23 - 26     Integer       resSeq       Residue sequence number
#  27          AChar         iCode        Code for insertion of residues
#  31 - 38     Real(8.3)     x            Orthogonal coordinates for X in A
#  39 - 46     Real(8.3)     y            Orthogonal coordinates for Y in A
#  47 - 54     Real(8.3)     z            Orthogonal coordinates for Z in A
#  55 - 60     Real(6.2)     occupancy    Occupancy
#  61 - 66     Real(6.2)     tempFactor   Temperature factor
#  77 - 78     LString(2)    element      Element symbol, right-justified
#  79 - 80     LString(2)    charge       Charge on the atom
#------------------------------------------------------------------------------
   pdb = open(file,"r")
   if (len(file) == 0):                                  # check for empty file
      print("ReadPDB: empty file ",file); exit(1)

   atoms = [Atom()]
   pdb.seek(0,0)
   for line in pdb:
      if (line[0:6] == "ATOM  " or line[0:6] == "HETATM"):
         atom = Atom()
         atom.name = line[12:16].strip()
         atom.resi = line[17:21].strip()
         atom.segm = line[21:22].strip()
         atom.ires = int(line[22:26])
         atom.r.x  = float(line[30:38])
         atom.r.y  = float(line[38:46])
         atom.r.z  = float(line[46:54])
         atom.occp = float(line[54:60]) if line[54:60] != "      " else 0e0
         atom.beta = float(line[60:66]) if line[60:66] != "      " else 0e0
         atom.symb = line[76:78].strip()
         atoms.append(atom)
   pdb.close()

   return atoms

#==============================================================================
def WritePDB(file, atoms):
#------------------------------------------------------------------------------
#  Writes a pdb file for atoms.
#------------------------------------------------------------------------------
   pdb = open(file,"w")
   
   strfrm = "ATOM  {:5d} {:4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}" + \
            "{:6.2f}{:6.2f}          {:2s}\n"
   natm = len(atoms) - 1
   for iatm in range(1,natm+1):
      atom = atoms[iatm]
      pdb.write(strfrm.format(iatm,atom.name,atom.resi,atom.segm[:1],
                              atom.ires,atom.r.x,atom.r.y,atom.r.z,
                              atom.occp,atom.beta,atom.symb))
   pdb.write("END\n")
   pdb.close()
