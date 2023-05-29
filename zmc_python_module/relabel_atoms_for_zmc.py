#!/opt/local/bin/python2.7

#############################
#
# Eric Chan september 2015
#
##------------------------------
#
# this is the original version 
# it relies on .cif format as the input file and will output a relabeled cif.
#
# the plan for v2.0 of this program is that it will create a centre of mass.
# in doing this it should be able to deal with multiple molecules which have been output 
# in a familar .mol2 format that is compatable with the original version of zmat_maker. 
#
# run by typing 
# $> python relabel_atoms_for_zmc.py rootname label
#
# where rootname is the rootname of the .cif 
# and label is the label of the atom you want to act as the center of rotation 
# for the Z-matrix in the ZMC program.
#
##################################################
 
import re
import sys
import os

import numpy as np
import openbabel

# import matplotlib.pyplot as plt
# import matplotlib

# import scipy
# import matplotlib  

############################################################
# read atoms from input cif and do the BFS to get the new atom labels 
############################################################

def readcif_and_bfs(infile,atomlabel_id_number):

 # now we use babel to get the data we want
 obConv = openbabel.OBConversion()
 obConv.SetInFormat("CIF")
 # obConversion.OpenInAndOutFiles("ethane.mol2","ethane.c

 mol = openbabel.OBMol()
 obConv.ReadFile(mol, infile)   # Open Babel will uncompress automatically

 # all you need to do is run BFS on the atom of intestest and the depth from the bonding atoms 
 # will tell you which partition it should be 

 # this is how to run the BFS algoritam 
 # read in the mol and the starting atom to return a tuple which contains the starting atom and 
 # the depth of the conectivity

 # important we should choose the starting atom based on what we want the final output to be  

 count = 0
 newid = []
 atomnums = [] 
 depths = [] 
 for bfs in openbabel.OBMolAtomBFSIter(mol,atomlabel_id_number): 
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

# print sorted_ids
# for i in range(np.size(newid)):
#  print "%5i %5i %5i %5i"%(newid[i], atomnums[i], depths[i], int(sorted_ids[i]))

##############################################################################

def get_atom_id_from_label(ciffile,label):

   atomlabel_id_number = 0
   with open(ciffile) as file:
        for line in file:
          line = line.replace("\n", " ") # get rid of newline 
          line = line.replace("\r", " ") # get rid of newline 
    # going down to the list of atoms
          if re.search('atom_site_label', line): break
        for line in file:
          if re.search('END', line): break
          if re.search('end', line): break

          line = line.replace("\n", "") # get rid of newline 
          line = line.replace("\r", "") # get rid of newline 
          match_underscore = re.search('^\_', line) 
          if not match_underscore :   # this is the block we want  
   #   now start counting 
   #        print line
           atomlabel_id_number += 1
           if re.search(label, line): break
                  
   return atomlabel_id_number 

################################################

from sys import argv
script, rootname, label = argv

ciffile = rootname+".cif"

atomlabel_id_number = get_atom_id_from_label(ciffile,label)

print "com atom id is "+str(atomlabel_id_number)

sorted_ids = readcif_and_bfs(ciffile,atomlabel_id_number)

# print ciffile
# ciffile = "186295B_relabel.cif"

####################################
# now we read in the cif 
#########################

cifheader = []
atomslist = []
atomlabel = []

with open(ciffile) as file:
        for line in file:
          line = line.replace("\n", " ") # get rid of newline 
          line = line.replace("\r", " ") # get rid of newline 

# goto the atoms in the list and save each line in atomlist[]
          cifheader.append(line)
          if re.search('atom_site_label', line): break
        for line in file:
          if re.search('END', line): break
          if re.search('end', line): break

          line = line.replace("\n", "") # get rid of newline 
          line = line.replace("\r", "") # get rid of newline 
          match_underscore = re.search('^\_', line) 
          if not match_underscore :   # this is the block we want 
           # print line
           atomslist.append(line)

# now just get the label and save in atomlabel[]
           asymbol = str(line.split()[1])   # in most .cif formats the symbol comes after the initial label 
         #  atomnum = re.sub('\D','',asymbol)   # \D match nondigits
           atomlabel.append(asymbol)
#           print atomnum 
          
# print atomslist
# print atomlabel

sorted_ids = np.array(sorted_ids)
atomslist = np.array(atomslist)
sortA = np.arange(60)+1.0       # this is just an array from 1-total num atoms for referencing 

# print sortA
print sorted_ids

# now read in the original cif and write to a new cif while relabeling on the fly 
# it should be possible to put in a dummy atom for the center of mass.

# this mean if your system has multiple molecules you have to cut the file into pieces and then relabel accordingly


outfile = rootname+"_relabel.cif"
filehandle = open(outfile, 'w')

# by using '\n' .join here we can print the entire array
filehandle.write('\n'.join(cifheader))
filehandle.write("\n") 
filehandle.write("_atom_site_type_symbol\n") 
filehandle.write("_atom_site_fract_x\n") 
filehandle.write("_atom_site_fract_y\n") 
filehandle.write("_atom_site_fract_z\n") 

# print ' '.join(atomslist[sorted_ids==i+1]) # tricky bit here with the ' '.join to get rid of the nparray ['  '] chars                
# here below we needed to split the string up and then reassemble it 

for i in range(np.size(atomslist)):
 line = ' '.join(atomslist[sorted_ids==i+1])
 words = line.split()
 label = words[0]
 new_label = re.sub('\d*','',label)
# print new_label+str(i+1)+' '+' '.join(words[1:])

 filehandle.write("%s \n"%(new_label+str(i+1)+' '+' '.join(words[1:])))

# filehandle.write("%s \n"%(' '.join(atomslist[sorted_ids==i+1]))) 

filehandle.write("#END \n") 
filehandle.close()

