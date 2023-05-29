#! /usr/bin/env python
#
# Support module generated by PAGE version 4.7
# In conjunction with Tcl version 8.6
#    Sep 09, 2016 03:11:38 PM
#
#  v1.0 - need to relpace all os.system commands with the subprocess module  
#    
#

import os
import sys
import re
import numpy as np
from tkFileDialog   import askopenfilename      
import ZMCGUI_functions as zf
import relabel_for_zmc_v2s1 
import subprocess 
from shutil import copyfile
from shutil import move

import manage_occupancies
import manage_contacts
import diffuse_input_controls
import ZMC_input_controls


try:
    from Tkinter import *
except ImportError:
    from tkinter import *

try:
    import ttk
    py3 = 0
except ImportError:
    import tkinter.ttk as ttk
    py3 = 1

def set_Tk_var():
    # These are Tk variables used passed to Tkinter and must be
    # defined before the widgets using them are created.
    global occfname
    global asize, bsize, csize, loc, zmat, nmol
    global latt_a, latt_b, latt_c, calpha, cbeta, cgamma
    global Sites

    occfname = StringVar()

    asize = StringVar()
    bsize = StringVar()
    csize = StringVar()
    loc = StringVar()
    zmat = StringVar()
    nmol = StringVar()
    Sites = StringVar()

    latt_a = StringVar()
    latt_b = StringVar() 
    latt_c = StringVar()
    calpha = StringVar()
    cbeta  = StringVar()
    cgamma = StringVar()

    asize.set(u"20")
    bsize.set(u"20")
    csize.set(u"20")
    loc.set(u"4")
    zmat.set(u"1")
    nmol.set(u"1")
    Sites.set(u"4")

    latt_a.set(u"10.0")
    latt_b.set(u"10.0")
    latt_c.set(u"10.0")
    calpha.set(u"90.0")
    cbeta.set(u"90.0") 
    cgamma.set(u"90.0")

    occfname.set(u"occ.txt")

    global ZMC_option_Diffuse
    ZMC_option_Diffuse = IntVar()
    ZMC_option_Diffuse.set(0)

    global DoRelabeling
    DoRelabeling = IntVar()
    DoRelabeling.set(0)

    global NoCentroid
    NoCentroid = IntVar()
    NoCentroid.set(0)


    global check_calcH
    check_calcH = IntVar()
    check_calcH.set(0)
    global check_calcC1
    check_calcC1 = IntVar()
    check_calcC1.set(0)


def BuildRandomOcc():
    print('ZMCGUI_support.BuildRandomOcc')
    sys.stdout.flush()
   
   # the purpose of newoccfname here is just to specify the path for creating the occupancy file 
    newoccfname = occfname.get() 
    if mol2file: newoccfname = re.sub('.mol2','_'+occfname.get(),str(mol2file))

    zf.make_random_occ(newoccfname,asize.get(),bsize.get(),csize.get(),loc.get(),zmat.get(),nmol.get())
    w.Text2.insert(END,"\n random occupancy file created with name "+occfname.get())
    w.Text2.see("end")

def ImportModelforZMC():
    print('ZMCGUI_support.ImportModelforZMC')
    sys.stdout.flush()
   # print str(mol2file)
    cparams = zf.get_cell_param_from_mol2(mol2file)
    
    latt_a.set(cparams[0])
    latt_b.set(cparams[1])
    latt_c.set(cparams[2])
    calpha.set(cparams[3])
    cbeta.set(cparams[4]) 
    cgamma.set(cparams[5])

  # decide to actually relabel the file or not but we will make a file with the tag _relabel anyway t 
    global mol2_relabel,cell_formula,non_hyd_atoms,atomtypes_total,num_hyd_atoms,total_atoms 
    mol2_relabel = re.sub('.mol2','_relabel.mol2',str(mol2file))
      
    if DoRelabeling.get() == 0 : 
         relabel_for_zmc_v2s1.relabel_for_zmc(str(mol2file))
         cell_formula = relabel_for_zmc_v2s1.get_cell_formula(str(mol2file))

    if DoRelabeling.get() == 1 : 
         copyfile(mol2file,mol2_relabel)
         cell_formula = relabel_for_zmc_v2s1.get_cell_formula(str(mol2file))


    w.Text2.insert(END,"\n cell formula : %s\n" %(cell_formula))
    w.Text2.see("end")

    print cell_formula
    formula_bad = re.search('\+',cell_formula)
    if formula_bad:
     w.Text2.insert(END," found a problem with formula  \n" )
     w.Text2.insert(END," numbers not properly assigned  \n" )
     w.Text2.see("end")
     atomtypes_total = 2     
     total_atoms = 16
     num_hyd_atoms = 8
     non_hyd_atoms = 8
    else:
     cfs = np.reshape(cell_formula.split(),(-1,2))   # separate str formula into columns 
     cfs_atomtypes = cfs[:,0]                        # separate atom types 
     cfs_natoms = cfs[:,1].astype(np.int)            # separate atom numbers 
     atomtypes_total =  np.size(cfs_atomtypes)     
     total_atoms = np.sum(cfs_natoms)
     num_hyd_atoms = cfs_natoms[cfs_atomtypes=='H']
     non_hyd_atoms = (total_atoms - num_hyd_atoms[0])
 
    w.Text2.insert(END,"\n total atoms: %i\n" %(total_atoms))
    w.Text2.insert(END,"\n atom types (all): %i\n" %(atomtypes_total))
    w.Text2.insert(END,"\n non-H atoms: %i\n" %(non_hyd_atoms))
    w.Text2.see("end")

    global projectname,qxyzfile,mapfile
    projectname = re.sub('.mol2','',str(mol2_relabel))
    qxyzfile = projectname+'.qxyz'
    mapfile = projectname+'.map'

    w.Text2.insert(END,"\n mol2 file sucessfuly converted")
    w.Text2.see("end")


def OpenDataFile():
    global mol2file,workpath,headername
    print('ZMCGUI_support.OpenDataFile')
    sys.stdout.flush()
    mol2file = askopenfilename() 
    print mol2file
 
 #  havent figured out how to get text box to overwirte 
 #  w.Text1.delete()

    w.Text1.insert(END,mol2file)
    w.Text1.see("end")

  # strips out the filename at the end of the directory tree  
  # this is sensitive to the operating system type i.e. linux or windows  
    mol2filename =  mol2file.split('/')[-1]
    workpath = re.sub(mol2filename,'',str(mol2file))
    headername = re.sub('.mol2','',mol2filename)

    w.Text2.insert(END,"\n current working directory : \n %s \n" %(workpath))
    w.Text2.insert(END,"\n header name : \n %s \n" %(headername))
    w.Text2.see("end")
    
def RunZMATMaker():
    print('ZMCGUI_support.RunZMATMaker')
    sys.stdout.flush()
    w.Text2.insert(END,"\n running: \n %s %s\n" %(zmatmaker_exe, mol2_relabel))
    zmatmaker_command_string = "%s %s" %(zmatmaker_exe, mol2_relabel)
    os.system(zmatmaker_command_string)

    num_sites, num_zmats, num_loc, num_type = zf.read_qxyz(qxyzfile) 

    Sites.set(num_sites)
    loc.set(num_loc)
    zmat.set(num_zmats)
    nmol.set(num_type)

    w.Text2.insert(END,"\n updated unit cell data\n")
    zf.make_map_file(mapfile,qxyzfile)
    w.Text2.insert(END,"\n map file created\n")
    w.Text2.see("end")

    if num_zmats > 1:
     for zmatn in range(num_zmats):
       zmatfile = projectname+'_'+str(zmatn+1)+'.zmat'
       if check_calcH.get() == 0: zf.edit_zmatrix(zmatfile,'H')
       if check_calcC1.get() == 0: zf.edit_zmatrix(zmatfile,'C1 ')
    else: 
     zmatfile = projectname+'.zmat'
     if check_calcH.get() == 0: zf.edit_zmatrix(zmatfile,'H')
     if check_calcC1.get() == 0: zf.edit_zmatrix(zmatfile,'C1 ')

def WriteZMCRunScript():
    print('ZMCGUI_support.WriteZMCRunScript')
    sys.stdout.flush()

    global ZMCinpf
    ZMCinpf = headername+"_relabel_ZMC.inp"  
    rsfh = open(workpath+'run_zmc_job.py','wb')

    myruntxt = """import os
import sys
import re
import numpy as np
import subprocess 
from shutil import copyfile
from shutil import move
"""
    rsfh.write(myruntxt)
    if ZMC_option_Diffuse.get() == 1: rsfh.write("zmccmd = \"%s --crystal --diffuse %s\"\n" %('ZMC',str(ZMCinpf))) 
    if ZMC_option_Diffuse.get() == 0: rsfh.write("zmccmd = \"%s --crystal %s\"\n" %('ZMC',str(ZMCinpf))) 
    rsfh.write("os.system(zmccmd)")
    rsfh.close()


def ExitZMCGUI():
    print('ZMCGUI_support.ExitZMCGUI')
    sys.stdout.flush()
    sys.exit() 

def GenQuickContacts():
    global contacts_outfile
    global contacts_trimmed,con2djgname
    global num_spring_types
    contacts_outfile = projectname+"_contacts_buck.txt" 
    contacts_trimmed = projectname+"_contacts_trimmed.txt"  
    con2djgname = projectname+"_contacts_fixed.all"  

    print('ZMCGUI_support.GenQuickContacts')
    sys.stdout.flush()
    w.Text2.insert(END,"\n generating contacts\n")
    gencontacts_command_string = "%s --min=0.0 --max=4.0 --num=10000 --collapse=0.0011 --collapsedistance --nosortkey %s > %s" %(contactsbuck_exe,mol2_relabel, contacts_outfile)
    os.system(gencontacts_command_string)

    con2djg_command_string = "%s %s %s > %s" %(con2djg_exe,mapfile,contacts_outfile,con2djgname)
    os.system(con2djg_command_string)
    w.Text2.insert(END,"\nSpring list convered for ZMC use\n")
    w.Text2.see("end")

  # gets rid of the negative valued springs and option for those to centriod  
    w.Text2.insert(END,"\n fixing spring list\n")
    w.Text2.see("end")
    num_spring_types = zf.fix_up_spring_list(projectname,contacts_outfile,con2djgname,contacts_trimmed,NoCentroid.get())  
    w.Text2.insert(END,"\nNumber of springs used : %i \n" %(num_spring_types))

    # num_spring_types = zf.read_contacts(contacts_trimmed)

    w.Text2.insert(END,"\n spring list done\n")
    w.Text2.see("end")

def GenQuickZMCinput():
    print('ZMCGUI_support.GenQuickZMCinput')
    sys.stdout.flush()

    global ZMC_inp_file
    ZMC_inp_file = projectname+"_ZMC.inp"  

    crysizepar =  [int(asize.get()),int(bsize.get()),int(csize.get()) ]
    cellpar =  [float(latt_a.get()),float(latt_b.get()),float(latt_c.get()),float(calpha.get()),float(cbeta.get()),float(cgamma.get())]

    w.Text2.insert(END,"\nGenerateing ZMC input file\n")

    zf.generate_ZMC_input_quick(projectname,headername,occfname.get(),ZMC_inp_file,crysizepar,cellpar,int(Sites.get()),int(zmat.get()),num_spring_types,contacts_trimmed)

    w.Text2.see("end")


def ManageOccupancies():
    print('ZMCGUI_support.ManageOccupancies')
    sys.stdout.flush()
   # called.create_Called(root, color="white", instance=1, geom= "+200+650")
    manage_occupancies.create_Manage_Occupancies(root, color="white", instance=1, geom= "+200+650")

def ManageContacts():
    print('ZMCGUI_support.ManageContacts')
    sys.stdout.flush()
    manage_contacts.create_Manage_Contacts(root, color="white", instance=1, geom= "+200+650")

def OpenDiffuseInputMenus():
    print('ZMCGUI_support.OpenDiffuseInputMenus')
    sys.stdout.flush()
    diffuse_input_controls.create_Diffuse_Input_Controls(root, color="white", instance=1, geom= "+200+650")

def RunZMC():
    print('ZMCGUI_support.RunZMC')
    sys.stdout.flush()
    print ZMC_option_Diffuse.get()

# edited so that it will work under linux
# note that subprocess arguments should not have any white space since this can cause errors 
    global ZMCinpf
    ZMCinpf = headername+"_relabel_ZMC.inp"  

    w.Text2.insert(END,"\n running ZMC\n")
    w.Text2.see("end")
   
    if ZMC_option_Diffuse.get() == 1:
     zmcrun = subprocess.Popen([zmc_exe,'--crystal','--diffuse', str(ZMCinpf) ],cwd=workpath)
     zmcrun.wait()

    if ZMC_option_Diffuse.get() == 0:
     zmcrun = subprocess.Popen([zmc_exe,'--crystal', str(ZMCinpf) ],cwd=workpath)
     zmcrun.wait()

    w.Text2.insert(END,"\n ZMC run completed\n")
    w.Text2.see("end")


def ZMCInputControl():
    print('ZMCGUI_support.ZMCInputControl')
    sys.stdout.flush()
    ZMC_input_controls.create_ZMC_Input_Controls(root, color="white", instance=1, geom= "+200+650")

def getvars():
    # this is a way we can pass variables from ZMC_GUI main interface  
    # into the other submodules
    # print 'in getvars'
    # print bin_dir,dif_input_gen_exe

    crysizepar =  [int(asize.get()),int(bsize.get()),int(csize.get()) ]
    cellpar =  [float(latt_a.get()),float(latt_b.get()),float(latt_c.get()),float(calpha.get()),float(cbeta.get()),float(cgamma.get())]


    zmc_gui_vars = []
    zmc_gui_vars.append(bin_dir)
    zmc_gui_vars.append(dif_input_gen_exe)
    zmc_gui_vars.append(workpath)  
    zmc_gui_vars.append(headername)    
    zmc_gui_vars.append(crysizepar)    
    zmc_gui_vars.append(cellpar)    
    zmc_gui_vars.append(total_atoms)
    zmc_gui_vars.append(atomtypes_total)
    zmc_gui_vars.append(non_hyd_atoms)

    return zmc_gui_vars

def getZMCINPvars():
    # passing variables for ZMC input menu and ZMC options

    crysizepar =  [int(asize.get()),int(bsize.get()),int(csize.get()) ]
    cellpar =  [float(latt_a.get()),float(latt_b.get()),float(latt_c.get()),float(calpha.get()),float(cbeta.get()),float(cgamma.get())]

 #   print projectname

    zmc_gui_vars = []
    zmc_gui_vars.append(bin_dir)
    zmc_gui_vars.append(dif_input_gen_exe)
    zmc_gui_vars.append(workpath)  
    zmc_gui_vars.append(headername)    
    zmc_gui_vars.append(crysizepar)    
    zmc_gui_vars.append(cellpar)    
    zmc_gui_vars.append(total_atoms)
    zmc_gui_vars.append(occfname.get())
    zmc_gui_vars.append(int(Sites.get()))
    zmc_gui_vars.append(int(loc.get()))
    zmc_gui_vars.append(int(zmat.get()))
    zmc_gui_vars.append(num_spring_types)
    zmc_gui_vars.append(contacts_trimmed)

    return zmc_gui_vars

def init(top, gui, *args, **kwargs):
    global w, top_level, root
    w = gui
    top_level = top
    root = top

    w.Text1.insert(END,"/path/to/.mol2 \n")
    w.Text2.insert(END,"\nZMC GUI output")
    w.Text2.insert(END,"\n---------------\n")

    # initialize paths to bin files 
    global GUI_dir, bin_dir
    GUI_dir = os.getcwd()             
    w.Text2.insert(END,"\n GUI path: \n %s \n" %(GUI_dir))

  #  bin_dir = os.path.dirname(os.getcwd()) + '\\bin\\'       
  # here we need to make an environ var to point to the executables  
    bin_dir = os.environ.get('ZMCGUI_BIN')         
    w.Text2.insert(END,"\n ZMC executables path: \n %s \n" %(bin_dir))

    global zmatmaker_exe, contactsbuck_exe,con2djg_exe,zmc_exe,dif_input_gen_exe
    zmatmaker_exe = bin_dir + 'zmat_maker'
    contactsbuck_exe = bin_dir + 'contacts_buck'
    con2djg_exe = bin_dir + 'con2djg'
    zmc_exe = bin_dir + 'ZMC'
    dif_input_gen_exe = bin_dir + 'diffuse_input_generator'
    print bin_dir

def destroy_window():
    # Function which closes the window.
    global top_level
    top_level.destroy()
    top_level = None

if __name__ == '__main__':
    import ZMCGUI
    ZMCGUI.vp_start_gui()

