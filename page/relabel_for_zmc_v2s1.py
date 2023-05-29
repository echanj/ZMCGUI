#!python

###################################
#
# This is subversion v2.1 (modified for use in the GUI)
# also added a routine to return cell formula  back into the GUI
# when reading in the .mol2 file  
#
# It is a modification to debug the automatic bond assignment problems associated with 
# when the centroid ends up being situated close to many other atoms which is often the case.
#
# currently set up for reading in a mol2 file that has all the atoms appropriotly packed.
# if starting from a .cif, read into mercury and then pack and then save as .mol2 
#
# run by typing 
# python relabel_for_zmc_v2.py  myStructure.mol2
# 
# this version tested and modified to be working on windows python running from cygwin terminal 
# there is a strange bug in  this version that does not allow the temporary .xyz file to be deleted within the script 
#
# read in a molecule using obabel and calcualte the center of mass (R) using
# R = 1/M*sum(m.r) where M is the sum of all masses. r should represent a vector component
# which means you have to calcualte for each vector component.
#
#  the centre of mass should satisfy the condition sum[m(r-R)] = 0
#
#####################################

import re
import sys
import os

import numpy as np
import openbabel


def get_atom_symbol(atom):

    atomicnum = atom.GetAtomicNum() # get the atomic number 

    if atom.IsHydrogen():
     atom_symbol = "H"  
    if atomicnum == 5:
     atom_symbol = "B"  
      
    if atom.IsCarbon():   
     atom_symbol = "C"  

    if atom.IsNitrogen():   
     atom_symbol = "N"  

    if atom.IsOxygen():  
     atom_symbol = "O"  

    if atomicnum == 9:
     atom_symbol = "F"  

    if atomicnum == 13:
     atom_symbol = "Al"  

    if atomicnum == 14:
     atom_symbol = "Si"  

    if atomicnum == 15:
     atom_symbol = "P"  

    if atomicnum == 16:
     atom_symbol = "S"  

    if atomicnum == 17:
     atom_symbol = "Cl"  
     FFtype = "Cl"

    if atomicnum == 31:
     atom_symbol = "Ga"  

    if atomicnum == 32:
     atom_symbol = "Ge"  

    if atomicnum == 33:
     atom_symbol = "As"  

    if atomicnum == 34:
     atom_symbol = "Se"  

    if atomicnum == 35:
     atom_symbol = "Br"  

    if atomicnum == 49:
     atom_symbol = "In"  

    if atomicnum == 50:
     atom_symbol = "Sn"  

    if atomicnum == 51:
     atom_symbol = "Sb"  

    if atomicnum == 52:
     atom_symbol = "Te"  

    if atomicnum == 53:
     atom_symbol = "I"  

    if atomicnum == 11:
     atom_symbol = "Na"  

    if atomicnum == 20:
     atom_symbol = "Ca"  

    if atomicnum == 22:
     atom_symbol = "Ti"  

    if atomicnum == 26:
     atom_symbol = "Fe"  

    if atomicnum == 30:
     atom_symbol = "Zn"  

    if atomicnum == 43:
     atom_symbol = "Tc"  

    if atomicnum == 44:
     atom_symbol = "Ru"  

#
# this is use of the "try" and exception error handling components of python 
# see http://stackoverflow.com/questions/1592565/determine-if-variable-is-defined-in-python 
#

    try:
      atom_symbol
    except NameError:
      print "there is no atomic symbol described for this atom using this routine"
      print "check for errors or modify script "
 #   else:
 #     print "atom symbol is defined."

    return atom_symbol


def readmol_and_bfs(mol,com_id_number):

 # important we should choose the starting atom based on what we want the final output to be  

 count = 0
 newid = []
 atomnums = [] 
 depths = [] 
 for bfs in openbabel.OBMolAtomBFSIter(mol,com_id_number): 
  atomindex = bfs[0].GetIdx()
  depth = bfs[1]
  count += 1

  newid.append(count) 
  atomnums.append(atomindex) 
  depths.append(depth)

 newid    = np.array(newid) 
 atomnums = np.array(atomnums) 
 depths   = np.array(depths) 

# print newid
# print atomnums
# print depths

# easiest way to sort one array against the other
# but we want to aviod using lists

 sort_atomnums=np.sort(atomnums)
 sorted_ids = np.zeros(np.shape(newid)) 
 for i in range(np.size(newid)):
  sorted_ids[i]=newid[atomnums==sort_atomnums[i]] 

 return sorted_ids 


def get_center_of_mass(mol):

 print mol.GetMolWt()

 sumx=0.0
 sumy=0.0
 sumz=0.0
 totalmass=0.0

# loop through all atoms 
 for obatom in openbabel.OBMolAtomIter(mol): 
   atomidx = obatom.GetIdx() 
   atom_symbol = obatom.GetId() 
   amass = obatom.GetAtomicMass() 
   x = obatom.GetX() 
   y = obatom.GetY() 
   z = obatom.GetZ() 
 
   sumx = sumx + x*amass 
   sumy = sumy + y*amass 
   sumz = sumz + z*amass

   totalmass = totalmass + amass 

#   print carts.GetX(), x , fracs.GetX()

#  mol.NumAtoms()

 com_x = sumx/totalmass 
 com_y = sumy/totalmass 
 com_z = sumz/totalmass 

# create a new atom for the center of mass 
 centroid = mol.NewAtom()
 centroid.SetAtomicNum(6)   # carbon atom
 centroid.SetVector(com_x, com_y, com_z) # coordinates
 carts = centroid.GetVector()

 return carts


def output_temp_xyz(mol,natoms,sorted_ids):

 fhout = open('temp_mol.xyz', 'w')
 fhout.write(str(natoms-1)+"\n")
 fhout.write("\n")
  
 for i in range(natoms):
  old_id = int(np.squeeze(np.where(sorted_ids==i+1))+1) 
  obatom = mol.GetAtom(old_id)       
  atomidx = obatom.GetIdx() 
  atom_symbol = get_atom_symbol(obatom)
  x = obatom.GetX()
  y = obatom.GetY()
  z = obatom.GetZ()
  if i !=0: fhout.write ("%s %.10f %.10f %.10f \n"%(atom_symbol,x,y,z))

 fhout.close()


def clean_up_mol2(fname,infile,natoms_list):

 outfile = re.sub('.mol2','_relabel.mol2',fname)
 fhout = open(outfile, 'w')
 nmols = np.size(natoms_list)

 with open(fname) as file:
       for line in file:
        if re.search('@<TRIPOS>CRYSIN',line):
         line = file.next() 
         line = line.replace("\n", " ") # get rid of newline 
         line = line.replace("\r", " ") # get rid of newline 
         cell_params_str = line
         print cell_params_str

 with open(infile) as file:
       for line in file:
        if re.search('@<TRIPOS>MOLECULE',line):
         fhout.write("@<TRIPOS>MOLECULE\n")
         fhout.write(file.next()) 
         line = file.next() 
         words = line.split()
         fhout.write("%s %s %i\n" %(words[0],words[1],nmols))
         fhout.write(file.next()) 
        elif re.search('GASTEIGER',line):
         fhout.write("GASTEIGER\n")
         fhout.write("*****\n")
         fhout.write("generated using ZMC relabeler v2.0\n")
        elif re.search('@<TRIPOS>SUBSTRUCTURE',line):
         fhout.write(line) 

         atomscount = 1  
         for s in range(np.size(natoms_list)):
          molnum = s+1
          fhout.write("    %i RES%i        %i GROUP             0 ****  ****    0\n" %(molnum,molnum,atomscount))  
          atomscount = atomscount+natoms_list[s]
         
        elif re.search('@<TRIPOS>CRYSIN',line):
         fhout.write(line) 
         fhout.write(cell_params_str+"\n") 
        else:
         fhout.write(line)

def fixmol2formerge(natoms,tempname,outname):
 
   print natoms
   fhout = open(outname, 'w')
  
   # get the relevant data, re-arange manually then write to the output file 
   mol2atomdata = []
   mol2bonddata = []
   with open(tempname) as file:
         for line in file:
          if re.search('@<TRIPOS>ATOM',line): break
         for i in range(natoms):
          line = file.next() 
          mol2atomdata.append(line) 
         for line in file:
          if re.search('@<TRIPOS>BOND',line): break
         for line in file:
          mol2bonddata.append(line) 

# now begin construction of output file

   numbonds = np.size(mol2bonddata)+1 

   fhout.write("@<TRIPOS>MOLECULE\n") 
   fhout.write("temp_mol.xyz\n") 
   fhout.write(" %i %i 0 0 0\n" %(natoms, numbonds)) 
   fhout.write("SMALL\n") 
   fhout.write("GASTEIGER\n\n") 
   # write the atoms 
   fhout.write("@<TRIPOS>ATOM\n")

 # print the last atom first  
   words = mol2atomdata[-1].split()
#   print '1 '+' '.join(words[1:6])+' 1 LIG1 0.0000'      # tricky - print syntax
   fhout.write('1 '+' '.join(words[1:6])+' 1 LIG1 0.0000\n')      # tricky - print syntax
   for i in range(natoms-1):
      words = mol2atomdata[i].split()
#      print str(i+2)+' '+' '.join(words[1:])      # tricky - print syntax
      fhout.write(str(i+2)+' '+' '.join(words[1:])+'\n')      # tricky - print syntax

   # write the bonds  an extra bond needs to be created
   fhout.write("@<TRIPOS>BOND\n") 
   
   fhout.write("     1    1    2    1\n") 
   for i in range(np.size(mol2bonddata)):
      words = mol2bonddata[i].split()
      bno = words[0]
      ba = words[1]
      bb = words[2]
      bt = words[3]
      fhout.write('     %i    %i    %i    %s\n' %(int(bno)+1, int(ba)+1, int(bb)+1, bt)) 

   fhout.close() 

def relabel_for_zmc(fname):

 print fname

 #  use babel to get the data we want
 obConv = openbabel.OBConversion()
 obConv.SetInFormat("mol2")

 mol = openbabel.OBMol()       
 obConv.ReadFile(mol, fname)   
 print mol.GetMolWt

 splitmol = mol.Separate()

 num_molecules = np.size(splitmol)

# print num_molecules

 natoms_list = [] 
 for n in range(num_molecules): 

  outname = 'mol'+str(n)+'.mol2'
  tempname = 'temp.mol2'
  mol1 = splitmol[n]

  carts = get_center_of_mass(mol1)
  print "center of mass in cartesisan coordinates"
  print "C0 C %.6f %.6f %.6f" %(carts.GetX(), carts.GetY(), carts.GetZ())
  com_XYZ=np.array([carts.GetX(), carts.GetY(), carts.GetZ()],"float64")
 
  com_id = mol1.NumAtoms()
  print "com atom id is "+str(com_id)
  sorted_ids = readmol_and_bfs(mol1,com_id)
 
  print sorted_ids

  natoms = mol1.NumAtoms()
  natoms_list.append(natoms) # keep a list of the number of atoms in each molecule
  output_temp_xyz(mol1,natoms,sorted_ids) # this outputs each molecule as a temp.xyz file  

# read the temp file back in 
  obConv.SetInFormat("xyz")
  mol1 = openbabel.OBMol()       
  obConv.ReadFile(mol1, 'temp_mol.xyz')   

# bug here   

# there is  a problem wiht how the automatic bond assigmnet variys when the newly input atoms 
# including the centriod actually affect the connectivity of the system 
# attempt to modify so that it outputs the centriod as a seperate object 

# its likely the bonding an order in the mol2 file needs to be fixed and all reordered  prior to subseqent merging the mol2 files 


## re-create the center of mass atom 
  centroid = mol1.NewAtom()
  centroid.SetId(1)   # make this the first atom  
  centroid.SetAtomicNum(6)   # carbon atom
  centroid.SetVector(com_XYZ[0], com_XYZ[1], com_XYZ[2]) # coordinates
#  mol1.AddBond(1, 2, 1)  
#  carts = centroid.GetVector()
  
# now send back out as mol2 format 
  obConv.SetOutFormat("mol2")
  obConv.WriteFile(mol1, tempname)

# it should be possible here to use an auxillary subroutine to fix up the .mol2 appropriotly prior to merging

  fixmol2formerge(natoms,tempname,outname)

# again problems cleaning up files in windows - not to serious
#  os.remove(tempname)

# diagnostic exit 
#  print " diagnositic exit"  
#  exit()

# we can just use the perlscript merge_mol2 just continuosly appending the molecules onto each other   
  if n == 0: # forgot to handle the case of only single molecule in unit cell
   merged1 = 'merged'+str(n)+'.mol2'
   cmdstr = 'merge_mol2_for_ZMC.pl mol0.mol2 > %s' %( merged1) 
   os.system(cmdstr)
  elif n == 1:
   merged1 = 'merged'+str(n)+'.mol2'
   cmdstr = 'merge_mol2_for_ZMC.pl mol0.mol2 %s > %s' %(outname, merged1) 
   os.system(cmdstr)
  else :
   merged1 = 'merged'+str(n-1)+'.mol2'
   merged2 = 'merged'+str(n)+'.mol2'
   cmdstr = 'merge_mol2_for_ZMC.pl %s %s > %s' %(merged1, outname, merged2) 
   os.system(cmdstr)

# delete unwanted files 
 obConv.CloseOutFile()
 for n in range(num_molecules): 
  outname = 'mol'+str(n)+'.mol2'
  os.remove(outname)

 for n in range(num_molecules-1): 
  merged = 'merged'+str(n)+'.mol2'
  os.remove(merged)

# now clean up the mol2 file so it works properly 

 print np.size(natoms_list)
 print natoms_list 

 infile = 'merged'+str(num_molecules-1)+'.mol2'
 print infile,fname 

 clean_up_mol2(fname,infile,natoms_list)
 os.remove(infile)


def get_cell_formula(fname):

 #  use babel to get the data we want
 obConv = openbabel.OBConversion()
 obConv.SetInFormat("mol2")

 mol = openbabel.OBMol()       
 obConv.ReadFile(mol, fname)   
 return mol.GetSpacedFormula()


if __name__ == "__main__":

 from sys import argv
 script, fname = argv

 relabel_for_zmc(fname)


 # 
 #
 # carts2 = get_center_of_mass(mol2)
 #
 # or loop through all atoms 
 # for obatom in openbabel.OBMolAtomIter(mol1): 
 #   atomidx = obatom.GetIdx() 
 #   atom_symbol = obatom.GetId() 
 #   amass = obatom.GetAtomicMass() 
 #   x = obatom.GetX()
 #   y = obatom.GetY()
 #   z = obatom.GetZ()
 #   resno = obatom.GetResidue()
 #   obatom.SetIdx(int(sorted_ids[atomidx-1]))
 #
 # to get around the sorting problem  
 # output temporary sorted  .xyz file in cart coords 
 # to do this you need to provide labels for the atoms so manually create lookup table for this  
 # and fill with most common atoms based on atomic number 
 #
 # reread back into this routine as .xyz and remove bonds from the centriod atom
 # connect to the second atom in list then   
 #
 # individual atoms in each mol2 also need to be labeled appropriotely
 #
 # there may be a problem in that if the centroid is very far away from the rest 
 # of the molecule when it is read back in it may come up as a seperate molecule 
 # so you need to put  a control in place for this 
 #
 # cat individul .mol2 molecules into the same .mol2 which can be operated by 
 # zmat_maker    
 #
 # get a single atom from the molecule by index (starts at 1)
 #  obatom = mol.GetAtom(i+1)       
 #  atomidx = obatom.GetIdx() 
 #  print atomidx,int(np.squeeze(np.where(sorted_ids==i+1))) 
 # 

 





