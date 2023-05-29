#! /usr/bin/env python

#
# set of GUI functions for ZMC 
# 

import os
import sys
import numpy as np
import re
from shutil import copyfile
from shutil import move

def make_random_occ(fname,csize1,csize2,csize3,loc,zocc,mocc): 

#
# convert the input to and integer 
#
 cella = int(csize1) 
 cellb = int(csize2)
 cellc = int(csize3) 
 loc = int(loc) 
 zocc = int(zocc) 
 mocc = int(mocc) 
#
# use np.random.random_integers() in order to get a number from 1-zocc
# fname = 'occ.txt'
#

 fh =  open(fname,'wb') 

 for a in range(cella):
  for b in range(cellb):
   for c in range(cellc):
    for l in range(loc):
     z = np.random.random_integers(1,zocc,1)
     m = np.random.random_integers(1,mocc,1)
     fh.write("%i %i %i %i %i %i\n"%(a+1, b+1, c+1, l+1, z[0], m[0]))

 fh.close() 


def get_cell_param_from_mol2(mol2name):

 with open(mol2name) as file:
        for line in file:
         if re.search('CRYSIN', line):
          line = file.next()  # this is the handy trick to jump the next line
          cellparams = line.split()
 
 return cellparams 


def load_contact_data(contacts_outfile,con2djgname):

# fileA = projectname+'_buck.txt' 
# fileB = projectname+'_fixed.all' 
 fileA = contacts_outfile
 fileB = con2djgname
 A = np.loadtxt(fileA) # the only colum we need from this is the sprconsts 
 B = np.loadtxt(fileB) # this is the con2djg fixed datafile
 M = np.c_[B,A[:,11]]   #  master contact data array

 return  M

def fix_up_spring_list(projectname,contacts_outfile,con2djgname,contacts_trimmed,NoCentroid):

 outname_sprcon = projectname+'_inp_sprcon.txt'
 outname_size = projectname+'_inp_sizef.txt'

 rawC = load_contact_data(contacts_outfile,con2djgname)  #  loading a master arrays 


#############################
# strip out all contacts 
# i.e. -ve or 0 valued sprint constants according to ericingham technique 
# connections to centroids
#
# the final contact array is labeled rawK
 rawK = rawC[rawC[:,13]>0.0]   #  select springs greater than 0.0

 if NoCentroid == 0:
  rawK = rawK[rawK[:,3]!=1]     #  select centriod as destination atom
  rawK = rawK[rawK[:,10]!=1]     #  select centriod as destination atom

###################
# This bit sorts out the renumbering of the contacts - taken from manage_contacts_for_ZMC.py 
##################
 ctypes = rawK[:,12]
 U,ind = np.unique(ctypes,return_inverse=True)
 ctypes = ind+1

 rawK[:,12] = ctypes 
###########################

 Crows,Ccol = np.shape(rawK) 
  
 print np.shape(rawK)

 fhout = open(contacts_trimmed,'wb')
 sprcon_out = open(outname_sprcon,'wb')
 sizef_out = open(outname_size,'wb')

 fhout.write('# ol   oz   om   oat   da   db   dc   dl   dz   dm  dat      dist     type\n')

 oldtype = 0
 for i in range(Crows):
  R = rawK[i,0:11].astype(int) # convert float to integer 
  dist = rawK[i,11]
  ctyp = int(rawK[i,12])
  sprcon = float(rawK[i,13])

  fhout.write('%4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i  %.8f  %4i\n' %(R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8],R[9],R[10],dist,ctyp))

  if oldtype != ctyp:
   oldtype = ctyp   
   sprcon_out.write('SPRCON  %.8f  %4i ! %.8f\n' %(sprcon,ctyp,dist))
   sizef_out.write('SIZE  1.0000  %4i ! %.8f\n' %(ctyp,dist))

 fhout.close()
 sprcon_out.close()
 sizef_out.close()

# move the trimmed list to the ZMC working file 
 move(contacts_trimmed,con2djgname)

 return int(np.max(rawK[:,12])) 

######################################

def make_map_file(mapfile,qxyzfile):
 
   mapfile_handle = open(mapfile, 'w') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes

   with open(qxyzfile) as myfile :
        for line in myfile :
         words = line.split()
         site = int(words[0]) 
         location =  int(words[1])
         zmat = int(words[2])
         ztype = int(words[3])
         mapfile_handle.write("%i %i %i %i \n" %(site, location, zmat, ztype)) 

  # changed from previous 
  # for i in range(num_residues):
  #  site,location,zmat,ztype = i+1,i+1,1,1 
  #  mapfile_handle.write("%i %i %i %i \n" %(site, location, zmat, ztype)) 
   
   mapfile_handle.close()

def read_qxyz(qxyzfile) :
     site_list = []
     loc_list = []
     zmat_list = []
     type_list = []
     with open(qxyzfile) as myfile :
          for line in myfile :
           linedata = line.split()
           site_list.append(int(linedata[0]))
           loc_list.append(int(linedata[1]))
           zmat_list.append(int(linedata[2]))
           type_list.append(int(linedata[3]))
     num_sites = max(site_list)           
     num_zmats = max(zmat_list)           
     num_loc   = max(loc_list)           
     num_type  = max(type_list)           
     return num_sites, num_zmats, num_loc, num_type 


def read_contacts(contacts_outfile) :
     type_list = []
     with open(contacts_outfile) as myfile :
          for line in myfile :
           match_header = re.search( '^\s*\#Num', line)   # this is how to get rid of the header 
           if not match_header:
            linedata = line.split()
            type_list.append(int(linedata[9]))
     num_spring_types = max(type_list)           
#     print max(type_list) 
     return int(num_spring_types)

def generate_ZMC_input_quick(projectname,headername,occfname,ZMC_inp_file,crysizepar,cellpar,num_residues,num_zmats,num_spring_types,contacts_trimmed):
 
 outname_sprcon = projectname+'_inp_sprcon.txt'
 print "\ntestline "+outname_sprcon
 con2djgname = headername+"_relabel_contacts_fixed.all"  
 occfname = headername+"_"+occfname
 ZMC_inp_file_handle = open(ZMC_inp_file, 'wb') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes
 spring_const = 1.0 
 dist = 1.0 
 cellpar = np.array(cellpar,dtype=np.float32)

# note that the block of text can be constructed without any spaces 
 ZMCinputtxt = """HEADER %s 
! Crystal geometry information
CRYSTAL %i %i %i
ZMATFILE  1  %s_relabel.zmat
QXYZFILE  1  %s_relabel.qxyz
CELL %.5f %.5f %.5f %.5f %.5f %.5f
CONTACTFILE %s
OCCFILE %s
!
! MC parameters
!
TEMPERATURE 1.0
MCCYCLES 3
XYZWIDTH 0.1
QWIDTH 0.10
INWIDTH 0.10
XYZINITW 0.0
QINITW 0.0
ININITW 0.0 
INCUPDATE 1
BADJUST 0 3.0
!
! For error checking, we can put some things in explicitely
!
NUMZMATS %d
NUMLOCS %d
NUMINSPRCON 0
NUMSPRCON %d
NUMINTERNAL ZMAT 1 0
NUMCROSS ZMAT 1 0

"""

 # print s %(spring_const,spring_const)
 ZMC_inp_file_handle.write(ZMCinputtxt %(headername,crysizepar[0],crysizepar[1],crysizepar[2],headername,headername,cellpar[0],cellpar[1],cellpar[2],cellpar[3],cellpar[4],cellpar[5],con2djgname,occfname,num_zmats,num_residues,num_spring_types))

# now just append the force constants file 
 with open(outname_sprcon) as myfile:
        for line in myfile:
          ZMC_inp_file_handle.write(line)

 ZMC_inp_file_handle.close()


# for the time being - in the full editing verison  
# we can generate the required ZMAT and QXYZ inputs by looing through automatically
# although we could provide user capability to edit explicitly 
# for now just setup means to meore easly change essential parameters so we can 
# work with both modulation and monte carlo using the GUI
# this means adding some switches to run ZMC 

def generate_ZMC_input_full(projectname,headername,occfname,ZMC_inp_file,crysizepar,cellpar,num_residues,num_zmats,num_spring_types,contacts_trimmed,ZMCinputVars,ZMC_modwave_options):
 
 runZMC_modwave,QModType,QVEC1,QVEC2,QVEC3,QPOL1,QPOL2,QPOL3,QZOCC,QCONC,QDIR1,QDIR2,QDIR3,QAMP = ZMC_modwave_options
 
 print runZMC_modwave
 print QVEC1
 
 outname_sprcon = projectname+'_inp_sprcon.txt'
 print "\ntestline "+outname_sprcon
 con2djgname = headername+"_relabel_contacts_fixed.all"  
 occfname = headername+"_"+occfname
 ZMC_inp_file_handle = open(ZMC_inp_file, 'wb') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes
 spring_const = 1.0 
 dist = 1.0 
 cellpar = np.array(cellpar,dtype=np.float32)

 ZMC_inp_file_handle.write("HEADER %s \n" %(headername))
 ZMC_inp_file_handle.write("! Crystal geometry information \n")
 ZMC_inp_file_handle.write("CRYSTAL %i %i %i \n" %(crysizepar[0],crysizepar[1],crysizepar[2]))
 
 if num_zmats > 1:
  for z in range(num_zmats):
    ZMC_inp_file_handle.write("ZMATFILE  %i  %s_relabel_%i.zmat \n" %(z+1,headername,z+1))
    ZMC_inp_file_handle.write("QXYZFILE  %i  %s_relabel.qxyz \n" %(z+1,headername))
 else:
    ZMC_inp_file_handle.write("ZMATFILE  1  %s_relabel.zmat \n" %(headername))
    ZMC_inp_file_handle.write("QXYZFILE  1  %s_relabel.qxyz \n" %(headername))

# note that the block of text can be constructed without any spaces 
 ZMCinputtxt = """CELL %.5f %.5f %.5f %.5f %.5f %.5f
CONTACTFILE %s
OCCFILE %s
!
! MC parameters
!
TEMPERATURE %.5f
MCCYCLES %i
XYZWIDTH %.5f
QWIDTH %.5f
INWIDTH %.5f
XYZINITW %.5f
QINITW %.5f
ININITW  %.5f
INCUPDATE %i
BADJUST %i %.5f
!
! For error checking, we can put some things in explicitely
!
NUMZMATS %d
NUMLOCS %d
NUMINSPRCON 0
NUMSPRCON %d
NUMINTERNAL ZMAT 1 0
NUMCROSS ZMAT 1 0

"""

 # print s %(spring_const,spring_const)
 ZMC_inp_file_handle.write(ZMCinputtxt %(cellpar[0],cellpar[1],cellpar[2],cellpar[3],cellpar[4],cellpar[5],con2djgname,occfname,ZMCinputVars[0],ZMCinputVars[1],ZMCinputVars[2],ZMCinputVars[3],ZMCinputVars[4],ZMCinputVars[5],ZMCinputVars[6],ZMCinputVars[7],ZMCinputVars[8],ZMCinputVars[9],ZMCinputVars[10],ZMCinputVars[11],ZMCinputVars[12],ZMCinputVars[13]))

 if runZMC_modwave == 1:  
  ZMC_inp_file_handle.write("SPRCON 1.000\n")
  ZMC_inp_file_handle.write("QMODTYPE %i \n" %(QModType))
  ZMC_inp_file_handle.write("QVECTOR %.5f %.5f %.5f \n" %(QVEC1,QVEC2,QVEC3))
  ZMC_inp_file_handle.write("QPOL %.5f %.5f %.5f  \n" %(QPOL1,QPOL2,QPOL3))
  ZMC_inp_file_handle.write("QAMP %.5f \n" %(QAMP)) 
  ZMC_inp_file_handle.write("QZOCC %i \n" %(QZOCC))
  ZMC_inp_file_handle.write("QCONC %.5f \n" %(QCONC))
  ZMC_inp_file_handle.write("QDIR %.5f %.5f %.5f  \n" %(QDIR1,QDIR2,QDIR3))
 
 else:

# now just append the force constants file 
  with open(outname_sprcon) as myfile:
        for line in myfile:
          ZMC_inp_file_handle.write(line)

  ZMC_inp_file_handle.close()




def edit_zmatrix(zmatfname,xatom):
     tempzmat = re.sub('.zmat','.zmattemp',zmatfname)
     tempzmatfh = open(tempzmat,'wb') 
     with open(zmatfname) as myfile:
        for line in myfile:
          tempzmatfh.write(re.sub(xatom,'x'+xatom,line))
     tempzmatfh.close()
     move(tempzmat,zmatfname)
# 
# contact_type_count = 0  # start counter at zero 
# 
# with open(contacts_trimmed) as myfile:
#       for line in myfile:
#        match_header = re.search( '^\s*\#Num', line)   # this is how to get rid of the header 
#
#        if not match_header:
#         line = line.replace("\n", " ") # get rid of newline 
#        # print line  
#        # print line , "\n" 
#         words = line.split()
#         dist = float(words[8]) 
#         contact_type = int(words[9]) 
#        # 
#        # spring_const = float(words[15]) # this was set to when sorkey was still output by the contacts_buck routine 
#         spring_const = float(words[11])
#         
#         if spring_const < 0: spring_const = 0.0  
#         # print spring_const
#
#         if contact_type != contact_type_count:  
#          contact_type_count += 1  # this is how to increment in python 
#          ZMC_inp_file_handle.write("SPRCON %6f  %3i !  %6f \n" %(spring_const, contact_type, dist))
#          print dist, contact_type, spring_const
#



if __name__ == '__main__':

 zmatfname =  '../testexample/aspirin/ACSALA07_relabel.zmat'
 edit_zmatrix(zmatfname,'H')

# print "select the function " 

# cparams =  get_cell_param_from_mol2("../ACSALA07.mol2")
# print cparams
# read_qxyz("../ACSALA07_relabel.qxyz")
#         
#         else:
#          line = file.next()
#          line = line.replace("\n", " ") # get rid of newline
#
#          note that the .replace() function is not as good as the re module
#          match_underscore = re.search('^\_', line)
#          newline = re.sub('^\w+',new_label,line)   # w is for word characters, W is for non-word chracters
#          cellparams = line.split()
#


