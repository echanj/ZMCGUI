#!python

#
# topas exporter 
# make sure the starting .cif is output from mercury 
# currently using v1.0 of the relabel for zmatmaker which does not have COM
# this restricts calculation to one molecule per cif file 
#
# run by typing 
#   python convert_for_topas.py  my_relabeled.cif my_relabeled.zmat my_relabeled.qxyz spacegroupsymbol
#
#
# eg.  python convert_for_topas.py  P44_HIVMI_relabel.cif P44_postDFT_cellfixed_riet.zmat P44_postDFT_cellfixed_riet.qxyz P21
# 
#
# need to get rid of potential for reading  of standard deviations when reading in .cif file  
#



import re
import sys
import os
import numpy as np

from sys import argv
script, cifin_name, zmatfile, qxyzfile, space_group_symbol = argv


def get_topas_header():

 TOPAS_header = """'' place header comment here 

xdd "c:\path\\to\data\here"

    One_on_X(@, 15116.32433)

	bkg @  88.1873152` -0.636830102`  1.06936095` -6.55227677`  4.17606429` -8.71863959` -8.00431967` -1.70097342`  0.271575302`  0.0993328722`  1.22066412` -2.11380948` -0.182308638` -2.99246236`  1.03584141` -2.38938205`  0.578013735`  0.937484999` -5.09592926`  0.86135728`  1.36333259`
	start_X  3
	finish_X  33
	LP_Factor( 90)
	Zero_Error(@, -0.01981`)
	Rp 217.5
	Rs 217.5

	lam
	   ymin_on_ymax  0.001
       la  0.653817 lo  1.540596 lh  0.501844
	   la  0.346183 lo  1.544493 lh  0.626579

'   Decompose(0.005)

'   Structure_Solution_Weighting


	str 
		CS_L(@, 9943.64796`_LIMIT_MIN_0.3)
	    TCHZ_Peak_Type(@, 1.99812`_LIMIT_MAX_2,@, -0.51928`,@,  0.04741`,, 0,@, 0.52729`,, 0)

	    ' MVW( 1327.987815, 1996.910684, 100)

		phase_name " place_structure_name_here "
		scale @  5.64835055e-005`

"""
 return(TOPAS_header)


def print_user_defined_options():

 user_defined = """\t\tUser_Defined_Dependence_Convolution(lor_fwhm, 1/Cos(Th), @, 0.01399)
\t\tUser_Defined_Dependence_Convolution(gauss_fwhm, 1/Cos(Th), @, 0.01514)
\t\tUser_Defined_Dependence_Convolution(hat, , @, 0.03890)
\t\tUser_Defined_Dependence_Convolution(circles_conv, -1/Tan(Th), @, 0.00045_LIMIT_MIN_0.0001)


"""
 return user_defined


# grab the parts of the cif file you need 

def get_cif_data(cifin_name,outfh):

 with open(cifin_name) as cif:
        for line in cif:
          line = line.replace("\n", " ") # get rid of newline 
          if re.search('_cell_length_a', line) :  outfh.write("\t\ta %s\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if re.search('_cell_length_b', line) : outfh.write("\t\tb %s\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if re.search('_cell_length_c', line) : outfh.write("\t\tc %s\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if re.search('_cell_angle_alpha', line) :  outfh.write("\t\tal %s\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if re.search('_cell_angle_beta', line) : outfh.write("\t\tbe %s\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if re.search('_cell_angle_gamma', line) : outfh.write("\t\tga %s\n\n\n"%( re.sub('\(\w*\)','',line.split()[1])) )
          if line.strip() == '_atom_site_label':  # Or whatever test is needed
           break
        atom_count = 0 
        for line in cif:
          line = line.replace("\n", " ") # get rid of newline 
          match_underscore = re.search('^\_', line) 
          if re.search("#END",line): break   
          if not match_underscore :   # this is the block we want 
            atom_count += 1 
            label = line.split()[0]
            symbol = line.split()[1]
            frac_x = line.split()[2] 
            frac_y = line.split()[3]
            frac_z = line.split()[4]
            frac_x = re.sub('\(\w*\)','',frac_x)
            frac_y = re.sub('\(\w*\)','',frac_y)
            frac_z = re.sub('\(\w*\)','',frac_z)
            if atom_count == 1 : # store the fractionals for the translation parameters 
             tx,ty,tz = frac_x,frac_y,frac_z
                       
            if symbol == 'H': 
             outfh.write("\t\tsite %s x %s y %s z %s occ %s 1 beq 1.0\n" %(label,frac_x,frac_y,frac_z,symbol) )               
            else:
             outfh.write("\t\tsite %s x %s y %s z %s occ %s 1 beq b%s 1.0\n" %(label,frac_x,frac_y,frac_z,symbol,symbol.lower()) )               

        print tx,ty,tz
        return tx,ty,tz


def read_write_zmatrix_info(zmatfile,outfh):

 outfh.write("\trigid\n") 

 with open(zmatfile) as zmat:
        zmat.next()  # skip the fisrt 2 lines
        zmat.next() 
        atom_labels = [] 
        atom_count = 0 
        for line in zmat:
         line = line.replace("\n", " ") # get rid of newline 
         label = line.split()[0]
         atom_2nd = line.split()[1]
         bond = line.split()[2]
         atom_3rd = line.split()[3]
         angle = line.split()[4]
         atom_4th = line.split()[5]
         torsion = line.split()[6]
         atom_labels.append(label)   # fill a labels list on the fly 
         atom_count += 1 
         if atom_count == 1 : 
          outfh.write("\t   z_matrix  %s \n" %(label))
         elif atom_count == 2 : 
          outfh.write("\t   z_matrix  %s \t %s \t RBD %s BDB(%.1f) \n" %(label, atom_labels[int(atom_2nd)-1], bond, float(bond)))
         elif atom_count == 3 : 
          outfh.write("\t   z_matrix  %s \t %s \t RBD %s BDB(%.1f) %s \t RAA %s TAB(%.1f) \n" %(label, atom_labels[int(atom_2nd)-1], bond, float(bond), atom_labels[int(atom_3rd)-1], angle, float(angle)))
         elif atom_count > 3 : 
          if float(torsion) > 179.0 and float(torsion) < 181.0 or float(torsion) > -1.0 and float(torsion) < 1.0 or float(torsion) < -179.0 and float(torsion) > -181.0  :
           outfh.write("\t   z_matrix  %s \t %s \t RBD %s BDB(%.1f) %s \t RAA %s TAB(%.1f) \t %s \t %s \n" %(label, atom_labels[int(atom_2nd)-1], bond, float(bond), atom_labels[int(atom_3rd)-1], angle, float(angle), atom_labels[int(atom_4th)-1], torsion))
          else:
           outfh.write("\t   z_matrix  %s \t %s \t RBD %s BDB(%.1f) %s \t RAA %s TAB(%.1f) \t %s \t RTA %s \n" %(label, atom_labels[int(atom_2nd)-1], bond, float(bond), atom_labels[int(atom_3rd)-1], angle, float(angle), atom_labels[int(atom_4th)-1], torsion))



 zmat_refinement_macros = """

\t   macro RBD { }
\t   macro BDB(v)	
\t\t   {min = v - .1; max = v + .1;}

\t   macro RAA {  }
\t   macro TAB(v)	
\t\t   {min = v - 10.0; max = v+10;}	  

\t   macro RTA { @ }
"""

 outfh.write(zmat_refinement_macros)


def write_rotate_translate(qxyzfile,outfh,tx,ty,tz):

 with open(qxyzfile) as qxyz:
        line = qxyz.readline()  # skip the line
        qr = float(line.split()[4])
        qi = float(line.split()[5])
        qj = float(line.split()[6])
        qk = float(line.split()[7])
        improper = line.split()[8]

# for formula for conversion of quaternion to euler angles 
# 
# see 
# https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
#
# these interpretations  still not thouroghly tested 
# also final qx qz qz values can be affected +/- 180 depending on the required input for topas 
#

 if improper == "F":
  Q = np.array([qr, qi, qj, qk])
 else:
  Q = np.array([qr, -qi, -qj, -qk])

 p1 =  2.0*(Q[0]*Q[1] + Q[2]*Q[3])
 p2 =  1.0 - 2.0*(Q[1]**2 + Q[2]**2) 
 t1 =  2.0*(Q[0]*Q[2] - Q[3]*Q[1])
 c1 =  2.0*(Q[0]*Q[3] + Q[1]*Q[2])
 c2 =  1.0 - 2.0*(Q[2]**2 + Q[3]**2) 

 phi = np.arctan2( p1,p2 )
 theta = np.arcsin( t1 )
 chi = np.arctan2( c1,c2 )

 qx = -phi*(180/np.pi)
 qy = -theta*(180/np.pi)
 qz = chi*(180/np.pi)+180

 oqx = -360+phi*(180/np.pi)
 oqy = theta*(180/np.pi)
 oqz = chi*(180/np.pi)


 outfh.write("\n\n")
 outfh.write("''if triclinic or monoclinic\n")
 outfh.write("\t   rotate @   %.6f qx 1\n"  %(qx) )
 outfh.write("\t   rotate @   %.6f qy 1\n"  %(qy))
 outfh.write("\t   rotate @   %.6f qz 1\n\n\n"  %(qz))

 outfh.write("\n\n")
 outfh.write("''if orthorhombic\n")
 outfh.write("\t   rotate @   %.6f qx 1\n"  %(oqx) )
 outfh.write("\t   rotate @   %.6f qy 1\n"  %(oqy))
 outfh.write("\t   rotate @   %.6f qz 1\n\n\n"  %(oqz))
 
 outfh.write("\t   translate\n" )
 outfh.write("\t      ta @  %.6f \n" %(float(tx)))
 outfh.write("\t      tb @  %.6f \n" %(float(ty)))
 outfh.write("\t      tc @  %.6f \n" %(float(tz)))



        
def main():

 rootname = cifin_name.replace(".cif", "")  

 topas_input_file = rootname+"_topas.inp"
 outfh = open(topas_input_file,'wb')
 
 outfh.write(get_topas_header())
 outfh.write("\n\n")
 outfh.write("\t\tspace_group  %s \n" %(space_group_symbol))

 tx,ty,tz = get_cif_data(cifin_name,outfh)
 outfh.write("\n\n\n")
 outfh.write(print_user_defined_options())

 # atomnum = 1;

 read_write_zmatrix_info(zmatfile,outfh)

 write_rotate_translate(qxyzfile,outfh,tx,ty,tz)

 outfh.write("\n\n\n")
 outfh.write(" normalize_FCs \n")
 outfh.write(" append_fractional \n")
 outfh.write(" append_bond_lengths \n") 
 outfh.write(" Out_CIF_STR(%s_rigid_zmat.CIF) \n" %(rootname))  
 outfh.write(" Out_FCF(%s_xdd.FCF) \n" %(rootname))
 outfh.write(" r_bragg  15.0000 \n")
 outfh.write(" view_structure \n")

if __name__ == "__main__":
  main()



###########################################################################################################
# spare code
#           line = line.replace("\n", " ") # get rid of newline 
#           match_underscore = re.search('^\_', line) 
#           if not match_underscore :   # this is the block we want 
#            words = line.split()
#            type_symbol = str(words[1])   # in most .cif formats the symbol comes after the initial label 
#            if type_symbol:
#             new_label = type_symbol+str(atomnum)+" "
#             atomnum += 1
#             newline = re.sub('^\w+',new_label,line)   # w is for word characters, W is for non-word chracters 
#             print newline               
#           else: print line
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

