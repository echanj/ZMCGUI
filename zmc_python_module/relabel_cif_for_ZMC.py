#!/usr/bin/python

# 
# The first atom that is read in will be considered atom "1" and all other atoms 
# will be assigned values according to the order they are read in from the file 
#
# if things dont work revert to a .cif format which is output by mercury program 
# make sure you remove the disorder from the file as well 
#
# The order of the atoms in the files needs to follow a dendritic path in order for zmatmaker to be of any use
# This is the original problem which is why I had sort for an atom labeling GUI extension independant of 
# other plugins simply for this reason. Development of the initial algorthams using programs such as materials studio     
# or jmol may be of some use or idealy this could be done by hand by clicking atoms sequentially.
#
# run by typing  send the output directly to the file you want 
#
# python  relabel_cif_for_ZMC.py in.cif > out.cif
#
# this script always giove an error warning which i have not figured out how to get around 
# it may be nessesary to use perl instead for regex, becaseu the python syntax can be frustrating
#
# there is no gernarl way to solve this problem using reg expression simple because most .cif formats are different 
#

import re
import sys
import os

from sys import argv
script, cifin_name = argv


print "run by typing  send the output directly to the file you want" 
print  "python  relabel_cif_for_ZMC.py in.cif > out.cif"
print  "you can ignore the errors but check to see things are working  "


# This is a standard way to loop through a block of data  
# using 'break'  
atomnum = 1;
with open(cifin_name) as cif:
        for line in cif:
          line = line.replace("\n", " ") # get rid of newline 
          print line               
          if line.strip() == '_atom_site_label':  # Or whatever test is needed
           break
        for line in cif:
          line = line.replace("\n", " ") # get rid of newline 
          match_underscore = re.search('^\_', line) 
          if not match_underscore :   # this is the block we want 
           words = line.split()
           type_symbol = str(words[1])   # in most .cif formats the symbol comes after the initial label 
           if type_symbol:
            new_label = type_symbol+str(atomnum)+" "
            atomnum += 1
            newline = re.sub('^\w+',new_label,line)   # w is for word characters, W is for non-word chracters 
            print newline               
          else: print line
        #  if not match_whitespace: print line
        #  else: break

        #  if line.strip() == '\#END':  
        #   break
        #  if line.strip() == 'loop_':  
        #   break
# print the rest of the cif out 
       # for line in cif:
       #   line = line.replace("\n", " ") # get rid of newline 
       #   print line               



