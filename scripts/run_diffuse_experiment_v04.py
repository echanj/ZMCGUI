
import re
import sys
import os
import numpy as np
# import square_ising 
import subprocess 
import glob

from shutil import copyfile, copy2
from shutil import move

# how to use copy file and move
## copyfile('temp.list', 'temp.txt')
# move('temp.list', contacts_trimmed)

####################################################
#
'''

v0.04 - (dec 2017) made minor workflow enhanements for DCDNB work

v0.03 - (oct 2016) minor changes to model ribbon system 

v0.02 - test use of python subprocess mdodule for handeling performing system calls in seperate directories
      - this might be a good way to spawn multiple diffuse runs
      - the usual way to do it is to use  os.chdir() and os.getcwd() but htis may be less appropriote for psudeo parallel operations 
      - http://stackoverflow.com/questions/13744473/command-line-execution-in-different-folder     
      - to use the glob module in the same manner as linux for wildcards.
        see http://stackoverflow.com/questions/11697245/copying-several-files-into-new-folder
      - https://jimmyg.org/blog/2009/working-with-python-subprocess.html

 apparently we need to be careful when using the program .wait() routine, but have not encountrerd any serious error so far.
    
 script to spawn diffuse experiment inputs from master files and run in seperate directory.
 the output pgm files will then be stored in a master pgm folder

 it should be possible to submit multiple experiments as seperate threads via an auxillary bash script 
 in this case it is experiment_design.sh  

'''
####################################################
#
# experiment parameters 
#

def setup_zmc_input(zmcinp,simsize,zmcP,occfile,exp_tag):

 asize,bsize,csize = simsize
 cycles,inspr,inwidth = zmcP

 fhout = open('temp.zmcinp'+str(exp_tag),'wb')
 
 with open(zmcinp) as myfile:
        for line in myfile:
         if re.match( 'CRYSTAL', line): 
          fhout.write('CRYSTAL %i %i %i\n' %(asize, bsize, csize)) 
         elif re.match( 'INSPR', line):
          words = line.split()
          fhout.write('INSPR '+str(inspr)+' '+' '.join(words[2:6])+'\n')
         elif re.match( 'INWIDTH', line):
          fhout.write('INWIDTH %.3f \n' %(inwidth)) 
         elif re.match( 'OCCFILE', line):
          fhout.write('OCCFILE %s \n' %(occfile)) 
         elif re.match( 'MCCYCLES', line):
          fhout.write('MCCYCLES %i \n' %(cycles)) 
 
         else: fhout.write(line)

 fhout.close()


def make_diffuse_inp(diffin,simsize,diffuse_params,exp_tag):

 asize,bsize,csize = simsize
 dest_diffinp =  re.sub('.in','_'+str(exp_tag)+'.in', diffin) 
 fhout = open(dest_diffinp,'wb')
 
 with open(diffin) as myfile:
        for line in myfile:
         # print line
         if re.search('simulation_size', line): 
          fhout.write('%i %i %i\n' %(asize, bsize, csize)) 
         elif re.search('Lot size', line): 
          fhout.write('%i %i %i\n' %(diffuse_params[0], diffuse_params[1], diffuse_params[2])) 
         elif re.search('Number of lots', line): 
          fhout.write('%i \n' %(diffuse_params[3])) 
         elif re.search('Number of atom sites per cell', line): 
          fhout.write('%i \n' %(diffuse_params[4])) 
         elif re.search('Number of atom types', line): 
          fhout.write('%i \n' %(diffuse_params[5])) 
         elif re.search('Subtract average lattice', line): 
          fhout.write('%s \n' %(diffuse_params[6])) 
         else: fhout.write(line)

 fhout.close()

 bin_name =  re.sub('.in','_'+str(exp_tag)+'.bin', diffin) 

# we no longer need to create input file because we can write this data directly to stdin using the subprocess module  

# fhout = open('input_diffuse_'+str(exp_tag),'wb')
# fhout.write('%s \n' %(dest_diffinp)) 
# fhout.write('%s \n' %(bin_name)) 
# fhout.close()
 
 return bin_name, dest_diffinp

def run_experiment(exp_tag,simsize,zmcP,paramJ,paramP,rootname,contacts_name,diffuse_calc_list,diffuse_params,pgm_mf,b2g_opt_list,b2g_opt):

# declare some general variables 

 zmcinp = rootname+'_ZMC.inp'
 diffuseout = rootname+'_ZMC.diffuse'
 
# create the new path 

 dir_path = './expfiles_'+str(exp_tag)+'/'

 if not os.path.exists(dir_path):
    os.makedirs(dir_path)

# copy over zmatricies and quaternions
 for zmat_file in glob.glob('*.zmat'):
    copy2(zmat_file,dir_path)
 for qxyz_file in glob.glob('*.qxyz'):
    copy2(qxyz_file,dir_path)

# copy over contacts file
 copy2(contacts_name,dir_path)

# run ising model -> pass variables using stdin  

 Ja,Jb,Jc,Jc2 = paramJ
 icycles,ianeal,iconc = paramP
# ilogfile = 'ising.log'
 ilogfile = open(dir_path+'ising.log', 'wb')
 
 isingrun = subprocess.Popen(['../occ_dcdnb' ],shell=False,cwd=dir_path,stdin=subprocess.PIPE,stdout=ilogfile)
 isingrun.stdin.write('%s \n' %(str(iconc)))
 isingrun.stdin.write('%s \n' %(str(Ja)))
 isingrun.stdin.write('%s \n' %(str(Jb)))
 isingrun.stdin.write('%s \n' %(str(Jc)))
 isingrun.stdin.write('%s \n' %(str(Jc2)))
 isingrun.stdin.write('%s \n' %(str(icycles)))
 isingrun.stdin.write('%s \n' %(str(ianeal)))
 isingrun.stdin.write('0 \n')
 isingrun.wait()
 ilogfile.close()
 
# setup other files and run programs
# occfile = 'occ_ising_'+str(exp_tag)+'.txt'
 occfile = 'occ_dcdnb.txt'
#  copy2(occfile,dir_path)
# move('ACSALA07_occ.txt',occfile)
# occfile = 'occ_ising.txt'

 setup_zmc_input(zmcinp,simsize,zmcP,occfile,exp_tag)

# make a backup of the occupancy file 
# bakoccfile = str(dir_path)+'occ_ising_'+str(exp_tag)+'.txt'
# move(occfile,bakoccfile)
# move(occfile,bakoccfile)

 dest_zmcinp = str(dir_path)+str(zmcinp)
 move('temp.zmcinp'+str(exp_tag), dest_zmcinp)
#  move('occ_ising.png','./png/occ_ising_'+str(exp_tag)+'.png')
#  move('corfunc.txt','./txt/corfunc_'+str(exp_tag)+'.txt')

# zmcrun = subprocess.Popen('ZMC --crystal --diffuse ',cwd=dir_path)

# edited so that it will work under linux
# note that subprocess arguments should not have any white space since this can cause errors 
#  exit()
 
 zmclogfile = open(dir_path+'zmc.log', 'wb')
 zmcrun = subprocess.Popen(['ZMC','--crystal','--diffuse', str(zmcinp) ],cwd=dir_path,stdout=zmclogfile)
 zmcrun.wait()
 zmclogfile.close()
 
# create input files for diffuse 

 dzmclogfile = open(dir_path+'diffuse.log', 'wb')
 b2glogfile = open(dir_path+'b2g.log', 'wb')
 for i in range(np.size(diffuse_calc_list)) :

  print diffuse_calc_list[i]  
  
  bin_name, dest_diffinp = make_diffuse_inp(diffuse_calc_list[i],simsize,diffuse_params,exp_tag)


# move diffuse input 
  move(dest_diffinp,str(dir_path)+dest_diffinp)

# move('input_diffuse_'+str(exp_tag),str(dir_path)+'input_diffuse_'+str(exp_tag))

# tricky way to pass required values to stdin
  dzmcrun = subprocess.Popen(['DZMC', str(diffuseout) ],shell=False,cwd=dir_path,stdin=subprocess.PIPE,stdout=dzmclogfile)
  dzmcrun.stdin.write('%s \n' %(dest_diffinp))
  dzmcrun.stdin.write('%s \n' %(bin_name))
  dzmcrun.wait()
 

#  os.system(dzmccmd)

  mystr = str(b2g_opt_list[b2g_opt[i]-1]).split()

#  print mystr
  if np.size(mystr) == 1 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],str(bin_name)],cwd=dir_path,stdout=b2glogfile)
  elif np.size(mystr) == 2 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],mystr[1],str(bin_name)],cwd=dir_path,stdout=b2glogfile)
  elif np.size(mystr) == 3 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],mystr[1],mystr[2],str(bin_name)],cwd=dir_path,stdout=b2glogfile)
  else :
   b2grun = subprocess.Popen(['bin2gray',str(bin_name)],cwd=dir_path,stdout=zmclogfile)

  b2grun.wait()
#  os.system(run_b2g)

# move the pgm to a master directory

  pgm_name =  re.sub('.bin','.pgm', bin_name) 
  move(dir_path+pgm_name,pgm_mf+pgm_name)

# keep the bin file in case you need to reprocess again
#  os.remove(bin_name)

 dzmclogfile.close()
 b2glogfile.close()

def main():

 from sys import argv
 script, exp_no = argv

 exp_no = int(exp_no)
# ncycles = 1 # int(ncycles)
# conc_A = float(conc_A)
 
# experiment arguments 

 simsize = [20,20,20]
 
 rootname = '300K_mixed_qeq_relabel'                        # transforming to N5  
 contacts_name = 'managed_contacts_fixed_revised.all'   # name of the contacts file
 pgm_mf = './pgm/'                            # pgm master folder - all pgm files from the experiemnt run will be moved here  

# diagnostic
 diffuse_params = [17,17,6,1,16,2,'n']       # list of input parameters for diffuse calc [lotsize1,lotsize2,lotsize3,nlots,natoms,atypes,subtract_ave_latt]             
# declare recipricol space projections
 diffuse_calc_list = []                       # create a list of diffuse scattering calcualtions based on master inputs that represent different scattering  
# diffuse_calc_list.append('diffuse_h0l.in')              
 diffuse_calc_list.append('diffuse_hk0.in')              
# diffuse_calc_list.append('diffuse_0kl.in')              
# diffuse_calc_list.append('diffuse_h1l.in')              
# diffuse_calc_list.append('diffuse_h2l.in')              
# diffuse_calc_list.append('diffuse_h3l.in')              
# diffuse_calc_list.append('diffuse_hk1.in')              
# diffuse_calc_list.append('diffuse_hk2.in')              
 diffuse_calc_list.append('diffuse_hk7.in')              

# declare bin2 gray options
 b2g_opt_list = []
 b2g_opt_list.append(' ')                                     # b2g option 1  
 b2g_opt_list.append(' --norm=20000 ')                        # b2g option 2 
 b2g_opt_list.append(' --hmirror  --vmirror  --norm=20000')                # b2g option 3  
 b2g_opt_list.append(' --twofold  --norm=20000 ')             # b2g option 4  
 b2g_opt = [3,3]    # int array  with size identical to diffuse_calc_list that contains the reference for which b2g option to use with which projection 
# b2g_opt = [3,3,3,3,3,3,3,3,3]    # int array  with size identical to diffuse_calc_list that contains the reference for which b2g option to use with which projection 
# b2g_opt = [4,3,3]    # int array  with size identical to diffuse_calc_list that contains the reference for which b2g option to use with which projection 

# zmc params
 zmcP = [1,10,0.001] # ncycles,insprcon,inwidth - the inwidth must be nonzero if any dihedral specified   

# note:  there is a tricky bug exists when using linux versions where if you specify that there are 
#        any dihedrals the program will segfault unless you put in a nonzero value for 
#        the inwidth 

# ising model params 
 paramJ =  [0.1,0.1,0.9,0.9]  # Ja,Jb,Jc,Jc2
# A =  [40,20] # asize,bsize 
 paramP =  [10,10,0.5] # ncycles,aneal,conc 

 run_experiment(exp_no,simsize,zmcP,paramJ,paramP,rootname,contacts_name,diffuse_calc_list,diffuse_params,pgm_mf,b2g_opt_list,b2g_opt)

if __name__ == "__main__":
    main()

