#!python

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
# Dec 2015 modified by Eric Chan to work better under the cygwin environment  
#
# Jan 2016 V2 - fix up organisation of spring list (get rid of negative valued springs as well and those coneting the centriod atoms)
#             - fix up creation of map file so that it matches the qxyz file 
# 
#######################################################################################################################
#
#
########################
# print cellpar  
#

import re
import sys
import os
import numpy as np 
from shutil import copyfile
from shutil import move

def generate_input_template(ZMC_inp_file,cellpar,num_residues,num_zmats,num_spring_types) :
 
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
 
 with open(contacts_trimmed) as myfile:
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
           ZMC_inp_file.write("SPRCON %6f  %3i !  %6f \n" %(spring_const, contact_type, dist))
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

def fix_up_spring_list():

# this fixes up the spring list so that only values with positive spring constant are included

 contact_type_count = 1  # start counter at zero 

 fhout = open('temp.list','wb') 

 with open(contacts_outfile) as myfile:
        line = myfile.readline()  # read in the first line of the file 
        line = myfile.readline()  # read in the second line of the file 
        olddist = float(line.split()[8]) 

 with open(contacts_outfile) as myfile:
        # read and write header
        line = myfile.readline() 
        fhout.write(line)  
        for line in myfile:
          words = line.split()
          oat = int(words[2]) 
          dat = int(words[7]) 
          dist = float(words[8]) 
          spring_const = float(words[11])
          if spring_const > 0.0 and not (dat == 1 or oat == 1) :           # get rid of the negative springs and those to centriod   
           if dist != olddist :  
            contact_type_count += 1  # this is how to increment in python 
            olddist=dist

           fhout.write("  "+'    '.join(words[0:9]) + "  %i  0.000  %0.4f\n" %(contact_type_count, spring_const))  

#          # print spring_const

 fhout.close()
# copyfile('temp.list', 'temp.txt')
 move('temp.list', contacts_trimmed)

#####
# here we begin running the program 
####

def main():

 global projectname
 global contacts_outfile
 global contacts_trimmed
 global mol2file 
 global qxyzfile
 global ZMC_inp
 global con2djgname 
 global mapfile

 from sys import argv
 # declare command line variables
 script, projectname = argv

 print "run by typing"
 print "python gen_ZMC_inputs.py projectname" 
   
 # the other way to do it is to use subprocess but it ends up being more complex 
 # for now just use os.system or exec() or execfile()
 # import subprocess
 
 # best thing to do is to declare inputs at the top

 # declare variables  
 # filename = sys.argv[1]

 contacts_outfile = projectname+"_contacts_buck.txt" 
 contacts_trimmed = projectname+"_contacts_trimmed.txt"  
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

 con2djg_command_string = "con2djg %s %s > %s" %(mapfile,contacts_trimmed,con2djgname)
 
 #####################
 # run zmatmaker and get contacts 
 #########################
 os.system(zmatmaker_command_string)
 os.system(gencontacts_command_string)
 ###################################### 

 fix_up_spring_list()  # this physicaly gets rid of the negative valued spring from the list prior to listing in the input file 

 ############################
 #  read the .mol2 and other files to get parameters for constructing the final input file  
 ############################ 

 cellpar,num_residues = read_mol2(mol2file)
 num_zmats = read_qxyz(qxyzfile)
 num_spring_types = read_contacts(contacts_trimmed)
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
 
 ZMC_inp_file = open(ZMC_inp, 'wb') # the w is for write mode, there is also 'a' for append and the '+' modifier..basically controls file attributes
 # filename = bulk_contacts_buck.txt 

 ################
 generate_input_template(ZMC_inp_file,cellpar,num_residues,num_zmats,num_spring_types)

 get_springs_from_list(ZMC_inp_file)

 ZMC_inp_file.close()

 #############

# you need to place this in order to run the function if this script is call 
# otherwise it just acts as a python module
if __name__ == "__main__":
  main()


