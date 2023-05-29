

import re
import sys
import os
import numpy as np
# import square_ising 


###############
# Author Eric J. Chan
# v1.0
# python scipt that can be used with pre existing diffuse input files to create the 
# nessesary input for smaller  recipricol space image boxes 
#
##############
# how to create the desired boxes 
#############
'''
1. designate a master difuse input file which has been created using the
   diffuse input generator that has an image resolution that is roughly similar to 
   what is desired.   
2. input the origin offset.  this is the bragg point for the origin (lower left hand of the image box) 
3. input desired rlu lengths along each axis of the box 
4.  adjust zoom, make large if small box    

example run 

python  edit_diffin.py diffuse_hk0.in

this will create diffuse_hk0_box.in 

''' 

from sys import argv
script, diffuse_inp = argv

######################
#
#        ^
# u-axis |
#        |
#  origin ----> v-axis
#
###############

# fhout = open(dest_diffinp,'wb')

# offset=[-4.05547,-4.05368,0.0] # origin h-offset,k-offset,l-offset to be expressed as rlu 
# rlu_rescale=[8.11094,8.10736]  # v-axis and u-axis in desired rlu


# 4 r.l.u. units of the upper left quadrant of the hk0 plane 
offset=[-4.0,0.0,0.0]  
rlu_rescale=[4.0,4.0]
zoom_scale=4.0 # rescale the number of pixels 

# offset=[-12.60476,0.0,-14.12445]  
# offset=[2.0,0.0,2.0]  
# rlu_rescale=[28.3815231425,27.95005]
# rlu_rescale=[4,3]
# zoom_scale=2.0 # rescale the number of pixels 
 
with open(diffuse_inp) as myfile:
        for line in myfile:
         if re.search('origin', line): 
           line = re.sub(',','',line)
           o_axis=np.array(line.split()[:3]).astype(np.float)
         #  print np.array(new).astype(np.float)
         #  new.astype(np.float)
         elif re.search('v_axis', line): 
           line = re.sub(',','',line)
           v_axis=np.array(line.split()[:3]).astype(np.float)
           v_pixels=np.array(line.split()[3]).astype(np.float)
         elif re.search('u_axis', line): 
           line = re.sub(',','',line)
           u_axis=np.array(line.split()[:3]).astype(np.float)
           u_pixels=np.array(line.split()[3]).astype(np.float)
         elif re.search('w_axis', line): 
           line = re.sub(',','',line)
           w_axis=np.array(line.split()[:3]).astype(np.float)
         elif re.search('lambda max', line): 
           stol_max=np.array(line.split()[0]).astype(np.float)
          # print line
           break

#       print o_axis 
#       print v_axis 
#       print u_axis 
#       print w_axis
#       print v_pixels,u_pixels,stol_max 
#       print "\n"

#       print (o_axis)-(o_axis)
#       print (v_axis)-(o_axis)
#       print (u_axis)-(o_axis)

rluv = np.linalg.norm((v_axis)-(o_axis))
rluu = np.linalg.norm((u_axis)-(o_axis))


rlupv = rluv/(v_pixels*zoom_scale) # rlu per pixel along v-axis
rlupu = rluu/(u_pixels*zoom_scale) # rlu per pixel along u-axis

pxav = np.rint(rlu_rescale[0]/rlupv) # pixels along v-axis
pxau = np.rint(rlu_rescale[1]/rlupu) # pixels along u-axis

re_rluv = np.rint(rlu_rescale[0]/rlupv)*rlupv  # rescaled v-axis in rlu
re_rluu = np.rint(rlu_rescale[1]/rlupu)*rlupu  # rescled u-axis in rlu
#
#       print (rluv) # reciprcol lattice units along v-axis 
#       print (rluu) # recipricol lattice units along u-axis
#       
#       print "\n"
#       print "copy and paste these coordinates into the diffuse input file\n"
#       print "\n"

outdiff = re.sub('.in','_box.in',diffuse_inp)

fhout = open(outdiff,'wb')

with open(diffuse_inp) as myfile:
        for line in myfile:
         if re.search('origin', line): break
         else: fhout.write(line) 
        line = myfile.next()  
        line = myfile.next()  
        line = myfile.next()  
        fhout.write( "%s %s ! origin \n" %(', '.join((((o_axis)-(o_axis))+offset).astype(np.str)) , str(', 0')))
        fhout.write( "%s %s ! v-axis \n" %(', '.join(((((v_axis-o_axis)/rluv)*re_rluv)+offset).astype(np.str)), ', '+str(int(pxav))))
        fhout.write( "%s %s ! u-axis \n" %(', '.join(((((u_axis-o_axis)/rluu)*re_rluu)+offset).astype(np.str)), ', '+str(int(pxau))))
        fhout.write( "%s %s ! w-axis \n" %(', '.join((((w_axis)-(w_axis))+offset).astype(np.str)), str(', 1')))
        for line in myfile:
         fhout.write(line) 

fhout.close()

#################################################

       #  line = file.next()  


  #       if re.search('simulation_size', line): 
  #        fhout.write('%i %i %i\n' %(asize, bsize, csize)) 
  #       elif re.search('Lot size', line): 
  #        fhout.write('%i %i %i\n' %(diffuse_params[0], diffuse_params[1], diffuse_params[2])) 
  #       elif re.search('Number of lots', line): 
  #        fhout.write('%i \n' %(diffuse_params[3])) 
  #       elif re.search('Number of atom sites per cell', line): 
  #        fhout.write('%i \n' %(diffuse_params[4])) 
  #       elif re.search('Number of atom types', line): 
  #        fhout.write('%i \n' %(diffuse_params[5])) 
  #       elif re.search('Subtract average lattice', line): 
  #        fhout.write('%s \n' %(diffuse_params[6])) 
  #       else: fhout.write(line)

# fhout.close()



