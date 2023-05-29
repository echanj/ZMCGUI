
import re
import sys
import os
import numpy as np
import ising3d 
import subprocess 
import glob

from shutil import copyfile, copy2
from shutil import move

# how to use copy file and move
## copyfile('temp.list', 'temp.txt')
# move('temp.list', contacts_trimmed)

####################################################
#
#
'''

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

def setup_zmc_input(zmcinp,asize,bsize,csize,cycles,inspr,inwidth,occfile,exp_tag):

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


def make_diffuse_inp(diffin,asize,bsize,csize,diffuse_params,exp_tag):

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

def run_experiment(exp_tag,asize,bsize,csize,Ja,Jb,Jc,Jcc,zmc_cycs,insprcon,inwidth,conc_A,rootname,contacts_name,diffuse_calc_list,diffuse_params,pgm_mf,b2g_opt_list,b2g_opt):

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

# get ising model parameters 
 J = [Ja,Jb,Jc,Jcc]

# these are ising model parameters that are currently fixed
 ncycles = 50
 aneal = 10

# run ising model 
 ising3d.ising(exp_tag,asize,bsize,csize,J,ncycles,aneal,conc_A)

# setup other files and run programs
 occfile = 'occ_ising_'+str(exp_tag)+'.txt'
 setup_zmc_input(zmcinp,asize,bsize,csize,zmc_cycs,insprcon,inwidth,occfile,exp_tag)

# make a backup of the occupancy file 
 bakoccfile = str(dir_path)+'occ_ising_'+str(exp_tag)+'.txt'
 move(occfile,bakoccfile)

 dest_zmcinp = str(dir_path)+str(zmcinp)
 move('temp.zmcinp'+str(exp_tag), dest_zmcinp)

# zmcrun = subprocess.Popen('ZMC --crystal --diffuse ',cwd=dir_path)

# edited so that it will work under linux
# note that subprocess arguments should not have any white space since this can cause errors 

 zmcrun = subprocess.Popen(['ZMC','--crystal','--diffuse', str(zmcinp) ],cwd=dir_path)
 zmcrun.wait()

# create input files for diffuse 

 dzmccmd = 'DZMC ' + str(diffuseout)  

 for i in range(np.size(diffuse_calc_list)) :

  print diffuse_calc_list[i]  
  
  bin_name, dest_diffinp = make_diffuse_inp(diffuse_calc_list[i],asize,bsize,csize,diffuse_params,exp_tag)

# move diffuse input 
  move(dest_diffinp,str(dir_path)+dest_diffinp)

# move('input_diffuse_'+str(exp_tag),str(dir_path)+'input_diffuse_'+str(exp_tag))

# tricky way to pass required values to stdin
  dzmcrun = subprocess.Popen(['DZMC', str(diffuseout) ],shell=False,cwd=dir_path,stdin=subprocess.PIPE)
  dzmcrun.stdin.write('%s \n' %(dest_diffinp))
  dzmcrun.stdin.write('%s \n' %(bin_name))
  dzmcrun.wait()

#  os.system(dzmccmd)

  mystr = str(b2g_opt_list[b2g_opt[i]-1]).split()

#  print mystr
  if np.size(mystr) == 1 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],str(bin_name)],cwd=dir_path)
  elif np.size(mystr) == 2 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],mystr[1],str(bin_name)],cwd=dir_path)
  elif np.size(mystr) == 3 :
   b2grun = subprocess.Popen(['bin2gray',mystr[0],mystr[1],mystr[2],str(bin_name)],cwd=dir_path)
  else :
   b2grun = subprocess.Popen(['bin2gray',str(bin_name)],cwd=dir_path)

  b2grun.wait()
#  os.system(run_b2g)

# move the pgm to a master directory

  pgm_name =  re.sub('.bin','.pgm', bin_name) 
  move(dir_path+pgm_name,pgm_mf+pgm_name)

# keep the bin file in case you need to reprocess again
#  os.remove(bin_name)

def main():

 from sys import argv
 script, exp_no, conc_A  = argv

 exp_no = int(exp_no)
 conc_A = float(conc_A)
 
# experiment arguments 
# run_experiment(exp_tag,asize,bsize,csize,Ja,Jb,Jc,Jcc,zmc_cycs,insprcon,inwidth,conc_A,rootname):
 
 rootname = 'shift_N5'                         # transforming to N5  
 contacts_name = 'merged_N5_F4_revised.all'    # name of the contacts file
 pgm_mf = './pgm/'                             # pgm master folder - all pgm files from the experiemnt run will be moved here  

# diagnostic
 diffuse_params = [5,5,2,1,864,4,'n']          # list of input parameters for diffuse calc [lotsize1,lotsize2,lotsize3,nlots,natoms,atypes,subtract_ave_latt]             
# diffuse_params = [13,11,2,192,864,4,'e']          # list of input parameters for diffuse calc [lotsize1,lotsize2,lotsize3,nlots,natoms,atypes,subtract_ave_latt]             

# declare recipricol space projections
 diffuse_calc_list = []                        # create a list of diffuse scattering calcualtions based on master inputs that represent different scattering  
 diffuse_calc_list.append('diffuse_h0l.in')              
 diffuse_calc_list.append('diffuse_0kl.in')              
 diffuse_calc_list.append('diffuse_hk0.in')              
 diffuse_calc_list.append('diffuse_h1l.in')
 diffuse_calc_list.append('diffuse_1kl.in')              

# declare bin2 gray options
 b2g_opt_list = []
 b2g_opt_list.append(' ')                                     # b2g option 1  
 b2g_opt_list.append(' --norm=20000 ')                        # b2g option 2 
 b2g_opt_list.append(' --hmirror  --vmirror ')                # b2g option 3  
 b2g_opt_list.append(' --twofold  --norm=20000 ')             # b2g option 4  
 b2g_opt = [3,3,3,3,3]                                 # int array  with size identical to diffuse_calc_list that contains the reference for which b2g option to use with which projection 

# diagnostic run 
 run_experiment(exp_no,8,8,5,1.0,1.0,0.0,0.0,10,20,1.0,conc_A,rootname,contacts_name,diffuse_calc_list,diffuse_params,pgm_mf,b2g_opt_list,b2g_opt)

# actual run 
# run_experiment(exp_no,38,38,38,1.0,1.0,0.0,0.0,2000,50,1.0,conc_A,rootname,contacts_name,diffuse_calc_list,diffuse_params,pgm_mf,b2g_opt_list,b2g_opt)

if __name__ == "__main__":
    main()

