#!/opt/local/bin/python

#######################################################################################################################
# this script generate most of the inputs required to run the ZMC program in order build a model for the crystal 
# and to drive the simulation 
#
# run by typing
#
# python gen_ZMC_inputs.py projectname 
#
# the program assumes you have no disorder in your system and will produce a basic set of working files 
#
# you need to have a valid mol2file avaliable that can run zmat_maker.exe
# the best way to test is to run zmat_maker.exe 
#
# in order to make things more prgamatic you should define the seperate parts of the program as functions 
# 
#######################################################################################################################
#
#
########################
# print cellpar  
print "run by typing"
print "python gen_ZMC_inputs.py projectname" 

import re
import sys
import os
from sys import argv
import numpy as np 

# declare command line variables
script, projectname = argv
 
def generate_input_template(ZMC_input_file,cellpar,num_residues,num_zmats,num_spring_types) :

 # filename = bulk_contacts_buck.txt 
 
 spring_const = 1.0 
 dist = 1.0 
 cellpar = np.array(cellpar,dtype=np.float32)

# note that the block of text can be constructed without any spaces 
 ZMCinputtxt = """HEADER %s 
! Crystal geometry information
CRYSTAL 20 20 20
ZMATFILE  1  %s.zmat
QXYZFILE  1  %s.qxyz
CELL %.5f %.5f %.5f %.5f %.5f %.5f
CONTACTFILE %s
OCCFILE occ.txt
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
 ZMC_inp_file.write(ZMCinputtxt %(projectname,projectname,projectname,cellpar[0],cellpar[1],cellpar[2],cellpar[3],cellpar[4],cellpar[5],con2djgname,num_zmats,num_residues,num_spring_types))

def get_springs_from_list(ZMC_inp_file) :
 
 contact_type_count = 0  # start counter at zero 
 
 with open(contacts_outfile) as myfile:
        for line in myfile:
         match_header = re.search( '^\s*\#Num', line)   # this is how to get rid of the header 

         if not match_header:
          line = line.replace("\n", " ") # get rid of newline 
         # print line  
         # print line , "\n" 
          words = line.split()
          dist = float(words[8]) 
          contact_type = int(words[9]) 
         # 
         # spring_const = float(words[15]) # this was set to when sorkey was still output by the contacts_buck routine 
          spring_const = float(words[11])
          
          if spring_const < 0: spring_const = 0.0  
          # print spring_const
 
          if contact_type != contact_type_count:  
           contact_type_count += 1  # this is how to increment in python 
           ZMC_inp_file.write("SPRCON %6f  ! %3i  %6f \n" %(spring_const, contact_type, dist))
           print dist, contact_type, spring_const

def read_mol2(mol2file) :
     with open(mol2file) as myfile :
          for line in myfile :
           if re.search('\s+\d', line):     # white space followed by digit 
            linedata = line.split()
            break
          for line in myfile :
           if re.search('CRYSIN', line):  
            break
          for line in myfile :
           cellparams = line.split()
   
     return cellparams,int(linedata[2])

def read_qxyz(qxyzfile) :
     zmat_list = []
     with open(qxyzfile) as myfile :
          for line in myfile :
           linedata = line.split()
           zmat_list.append(int(linedata[2]))
     num_zmats = max(zmat_list)           
     return int(num_zmats)

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

def make_map_file(mapfile,num_residues):
 
   mapfile_handle = open(mapfile, 'w') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes

   for i in range(num_residues):
    site,location,zmat,ztype = i+1,i+1,1,1 
    mapfile_handle.write("%i %i %i %i \n" %(site, location, zmat, ztype)) 
   
   mapfile_handle.close()

#####
# here we begin running the program 
####

# you need to place this in order to run the function if this script is call 
# otherwise it just acts as a python module
if __name__ == "__main__":
   
 # the other way to do it is to use subprocess but it ends up being more complex 
 # for now just use os.system or exec() or execfile()
 # import subprocess
 
 # best thing to do is to declare inputs at the top

 # declare variables  
 # filename = sys.argv[1]
 contacts_outfile = projectname+"_contacts_buck.txt" 
 mol2file = projectname+".mol2" 
 qxyzfile = projectname+".qxyz" 
 ZMC_inp = projectname+"_ZMC.inp"
 con2djgname = projectname+"_contacts_fixed.all"
 mapfile=projectname+".map"
 
 print contacts_outfile
 print mol2file

 # os.system("contacts.exe --min=0.0 --max=4.0 --num=10000 --collapse=0.0011 --collapsedistance urea_merc.mol2  > bulk_contacts.txt")
 # subprocess.call(["contacts_buck.exe", "urea_merc.mol2"])
 # subprocess.Popen('echo hello world > new.txt',shell=True)
 
 zmatmaker_command_string = "zmat_maker %s" %(mol2file)
 gencontacts_command_string = "contacts_buck --min=0.0 --max=4.0 --num=10000 --collapse=0.0011 --collapsedistance --nosortkey %s > %s" %(mol2file, contacts_outfile)

 con2djg_command_string = "con2djg %s %s > %s" %(mapfile,contacts_outfile,con2djgname)
 
 #####################
 # run zmatmaker and get contacts 
 #########################
 os.system(zmatmaker_command_string)
 os.system(gencontacts_command_string)
 ###################################### 

 ############################
 #  read the .mol2 and other files to get parameters for constructing the final input file  
 ############################ 

 cellpar,num_residues = read_mol2(mol2file)
 num_zmats = read_qxyz(qxyzfile)
 num_spring_types = read_contacts(contacts_outfile)
# print cellpar, num_residues 

#############################################
# at this stage we should create a standard .map file
# since we know the how many molecules are in the unit cell   
###########################################

 make_map_file(mapfile,num_residues)

 # run con2djg 
 os.system(con2djg_command_string)

 #############################################
 # begin the part to generate the input file 
 #############################################
 
 ZMC_inp_file = open(ZMC_inp, 'w') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes
 
#################
 generate_input_template(ZMC_inp_file,cellpar,num_residues,num_zmats,num_spring_types)
 get_springs_from_list(ZMC_inp_file)

 ZMC_inp_file.close()
##############


