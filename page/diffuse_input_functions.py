#! /usr/bin/env python

#
# set of GUI functions for diffuse 
# 

import os
import sys
import numpy as np
import re
from shutil import copyfile
from shutil import move

def generate_diffuse_input(difgen_var,workpath,DIG_fname,cellpar,crysizepar):

 subtract_bragg, lot_a, lot_b, lot_c, nlots,nsites, ntypes, check_doPBC, checkCarbon, checkOxygen,checkNitrogen, checkSulfur, checkChlorine, checkFluorine, PixelsX, PixelsY, PixelsZ = difgen_var
 
# read in the previous diffuse file
 with open(workpath+DIG_fname) as myfile:
          content = myfile.readlines()


 scatf_C="""'C '      
2.31,20.8439,1.02,10.2075
1.58860,0.5687,0.865,51.6512,0.2156
0.0,0.0                   !fprime, f-double-prime
"""

 scatf_O="""'O '                 
3.0485,13.2771,2.2868,5.7011
1.5463,0.3239,0.867,32.9089,0.2508
0.0,0.0                   !fprime, f-double-prime
"""

 scatf_N="""'N'
12.2126,0.00570,3.13220,9.89330
2.01250,28.9975,1.16630,0.58260,-11.529
0.0,0.0
"""

 scatf_H="""'H '                           !
0.0000,0.0000,0.0000,00.0000
0.0000,00.0000,0.0000,0.0000,000.000
0.0,0.0                   !fprime, f-double-prime
"""

 scatf_S="""'S'
6.9053,1.4679,5.2034,22.2151
1.4379,0.2536,1.5863,56.172,0.8669
0.0,0.0
"""

 scatf_Cl="""'Cl'
11.4604,0.0104,7.1964,1.1662
6.2556,18.5194,1.6455,47.7784,-9.5574
0.0,0.0
"""

 scatf_F="""'F'
3.5392,10.2825,2.6412,4.2944
1.517,0.2615,1.0243,26.1476,0.2776
0.0,0.0
"""

 if check_doPBC == 0 : pbcstr = 'n'
 if check_doPBC == 1 : pbcstr = 'y'

 dfofh = open(workpath+'temp.diffuse.in','wb')
 dfofh.write(content[0])
 dfofh.write(content[1])
 dfofh.write(content[2])
 dfofh.write("%s %s %s  ! simulation size \n" %(crysizepar[0],crysizepar[1],crysizepar[2]))
 dfofh.write("%s \n"%(pbcstr))
 dfofh.write(content[5])

# had to chop this up to be able to adjust pixels 
# dfofh.write(content[6])
 v_axis =  content[6].split()
 dfofh.write(' '.join(v_axis[0:3])+'  '+str(PixelsX)+'  '+' '.join(v_axis[4:6])+'\n')
# dfofh.write(content[7])
 u_axis =  content[7].split()
 dfofh.write(' '.join(u_axis[0:3])+'  '+str(PixelsY)+'  '+' '.join(u_axis[4:6])+'\n')
# dfofh.write(content[8])
 w_axis =  content[8].split()
 dfofh.write(' '.join(w_axis[0:3])+'  '+str(PixelsZ)+'  '+' '.join(w_axis[4:6])+'\n')

 dfofh.write(content[9])
 dfofh.write("%i,%i,%i  ! Lot size  \n" %(lot_a,lot_b,lot_c)) 
 dfofh.write("%i     ! Number of lots  \n" %(nlots))
 dfofh.write("%i   ! Number of atom sites per cell  \n" %(nsites))
 dfofh.write("%i   ! Number of atom types  \n" %(ntypes))
 dfofh.write("%s   ! Subtract average lattice?  \n" %(subtract_bragg))
 if checkCarbon==1: dfofh.write(scatf_C)  
 if checkOxygen==1: dfofh.write(scatf_O)  
 if checkNitrogen==1: dfofh.write(scatf_N)  
 if checkSulfur==1: dfofh.write(scatf_S)  
 if checkChlorine==1: dfofh.write(scatf_Cl)  
 if checkFluorine==1: dfofh.write(scatf_F)  
 dfofh.write(scatf_H)
 dfofh.close()
    
 move(workpath+'temp.diffuse.in',workpath+DIG_fname)

def read_pgm(filename, byteorder='>'):
    """Return image data from a raw PGM file as numpy array.

    Format specification: http://netpbm.sourceforge.net/doc/pgm.html

    """
    with open(filename, 'rb') as f:
        buffer = f.read()
    try:
        header, width, height, maxval = re.search(
            b"(^P5\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n]\s)*)", buffer).groups()
    except AttributeError:
        raise ValueError("Not a raw PGM file: '%s'" % filename)
    # dt = numpy.dtype(numpy.uint16)
    dt = np.dtype(np.uint16)
    dt = dt.newbyteorder('>')
    return np.frombuffer(buffer,
                            # dtype='u1' if int(maxval) < 256 else byteorder+'u2',
                            dtype='u1' if int(maxval) < 256 else dt,
                            # dtype=dt,
                            count=int(width)*int(height),
                            offset=len(header)
                            ).reshape((int(height), int(width)))



if __name__ == '__main__':

 print "select the function " 
 workpath = '../testexample/aspirin/'
 DIG_fname = 'diffuse_h0l.inp'

 generate_diffuse_input('difgen_var',workpath,DIG_fname,'cellpar','crysizepar')

# cparams =  get_cell_param_from_mol2("../ACSALA07.mol2")
# print cparams
# read_qxyz("../ACSALA07_relabel.qxyz")



