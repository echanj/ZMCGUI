
# square lattice ising spin model (only 2-states)

# this is a standard square ising model which is
# based on fliping the sign of the spin. 

# we also calculate the spin correlation function 


# just based on what they are saying we will only consider the cell consiting of two chains 
# and only use the ising model to got occuplanies for the first chain  

import numpy as np
import pylab
import time

# subroutine to write to occupancy file 
def occupancy(A,asize,bsize,cella,cellb,cellc):

# adopt a strategy to ocupations on the square grid to specific locations 

# this is a tricky part of the procedure

# loc = 8
# m=1

 occ_fh = open("occ_ising.txt","w")

 for a in range(cella):
  for b in range(cellb):
   for c in range(cellc):

       z = 1 
    # mapping to locations here
       c1 = c*2
       c2 = c*2+1 
     #  standard format    
     #  occ_fh.write("%3i %3i %3i %3i %3i %3i \n" %(a+1, b+1, c+1, l, z, m))

       occ_fh.write("%3i %3i %3i %3i %3i %3i \n" %(a+1, b+1, c+1, 1, z, ((A[a,c1]-1)/2)+2))
       occ_fh.write("%3i %3i %3i %3i %3i %3i \n" %(a+1, b+1, c+1, 2, z, ((A[a,c1]-1)/2)+2))
       occ_fh.write("%3i %3i %3i %3i %3i %3i \n" %(a+1, b+1, c+1, 3, z, ((A[a,c2]-1)/2)+2))
       occ_fh.write("%3i %3i %3i %3i %3i %3i \n" %(a+1, b+1, c+1, 4, z, ((A[a,c2]-1)/2)+2))

 occ_fh.close()


# energy subroutine
def energy(J,A,xx,yy,bnd1,bnd2):

   Ja,Ja2,Jb = J

   px = bnd1[xx+1] 
   mx = bnd1[xx-1] 
   px2 = bnd1[xx+2] 
   mx2 = bnd1[xx-2] 
   py = bnd2[yy+1] 
   my = bnd2[yy-1]

   e =  -Ja*A[xx,yy]*A[px,yy] 
   e = e-Ja*A[xx,yy]*A[mx,yy] 
   e = e-Ja2*A[xx,yy]*A[px2,yy] 
   e = e-Ja2*A[xx,yy]*A[mx2,yy] 
   e = e-Jb*A[xx,yy]*A[xx,py] 
   e = e-Jb*A[xx,yy]*A[xx,my] 

   return e

# define a spin correlation function here 

def correlation(A,asize,bsize,c1,c2):

 bnd1 = np.array([np.arange(asize),np.arange(asize)])
 bnd2 = np.array([np.arange(bsize),np.arange(bsize)])
 bnd1 = bnd1.reshape(asize*2)
 bnd2 = bnd2.reshape(bsize*2)

 C = np.zeros((c1,c2)) # declare correlation function range 

 # get the average spin 
 # the mean value for the spin 
 sigma0 = np.mean(A)  

 #  px = bnd1[xx+1] 
 #  mx = bnd1[xx-1] 
 #  py = bnd2[yy+1] 
 #  my = bnd2[yy-1]

 for ca in range(c1):
  for cb in range(c2):
   sum_px = 0.0 
   for a in range(asize):
    for b in range(bsize):
     px = bnd1[a+ca] 
     py = bnd2[b+cb] 
     sum_px =  sum_px + A[a,b]*A[px,py]
   C[ca,cb] = (sum_px/(asize*bsize))-(sigma0**2) 

 sum_px = 0.0 
 sum_pxx = 0.0 
 for a in range(asize):
  for b in range(bsize):
   px = bnd1[a+1] 
   pxx = bnd1[a+2] 
   sum_px =  sum_px + A[a,b]*A[px,b]
   sum_pxx =  sum_pxx + A[a,b]*A[pxx,b]

 corr = (sum_px/(asize*bsize))-(sigma0**2) 
 corr2 =  (sum_pxx/(asize*bsize))-(sigma0**2) 
 print "corr = %.6f" %(corr)
 print "corr2 = %.6f" %(corr2)
 print "ave sigma0 = %.6f" %(np.mean(A))
 print "ave sigma0^2 = %.6f" %(sigma0**2)
 print sum(A)
 
# this is how to adjust np.array print options 
 np.set_printoptions(precision=3)
# np.set_printoptions(threshold=np.nan)
# np.set_printoptions(edgeitems=2)
 np.set_printoptions(linewidth=300)

# print out the correlation function 
 print C
 mycorr =  open('corfunc.txt','w')
 for ca in range(c1):
  for cb in range(c2):
   mycorr.write(" %6f" %(C[ca,cb]))
  mycorr.write("\n")


# now from the correlation we can define the probability matrix that a certain spin will occur
# this is based on eq 2.15 in richards book. 
  
 P = np.zeros((c1,c2)) # declare correlation function range 
 ma = np.sum(A==1.0)/float(np.size(A))  # get the proportion of spin up     
 mb = np.sum(A==-1.0)/float(np.size(A))  # get the proportion of spin down 
 
 for pa in range(c1):
  for pb in range(c2):
   # P[pa,pb] = C[pa,pb] 
    P[pa,pb] = ma**2+(C[pa,pb]*(ma*mb))
 
 print "proportion spin up = %.6f " %(ma)
 print "proportion spin down = %.6f " %(mb)
 mycorr.write("proportion spin up = %.6f " %(ma))
 mycorr.write("proportion spin down = %.6f " %(mb))

 print P 
 print np.size(A)
 mycorr.close()

# return C

###############

def ising(L,J,P):

 asize,bsize = L
 ncycles,aneal,temp = P
 
# this is how to create a random 2D spin magnet array all in one line
 A = np.random.random_integers(0,1,(asize,bsize))*2-1


# but since we want exactly 50% occupancy we need to constrain the array 
# create 10 by 5 arrays then concatenate 
# note that the concatenate function is sensitive to the axis choice
# a = -np.ones((asize,bsize/2))
# b = np.ones((asize,bsize/2))
# A =np.concatenate((a,b),axis=1)

# now in order to randomly shuffle the array in python you need to make it one dimensional  
# A=A.reshape(asize*bsize)
# np.random.shuffle(A)
# np.random.shuffle(A)
# A=A.reshape(asize,bsize)

# ########################################################
# # extra bit here fopr ploting 
# ##############################################
# 
# # the sum of the spins should always be zero
# # you can plot the result in black and white using 
# pylab.plt.imshow(A,cmap='Greys',interpolation='nearest')
# pylab.plt.ion() # this turns ploting interactive mode on so we can update frames  
# pylab.plt.show()
# # reshuffle the array and redraw 
# 
# time.sleep(0.05)
# for i in range(10): 
#  A=A.reshape(asize*bsize)
#  np.random.shuffle(A)
#  A=A.reshape(asize,bsize)
#  pylab.plt.imshow(A,cmap='Greys',interpolation='nearest')
#  pylab.plt.draw()
#  time.sleep(0.10)
# 
# pylab.plt.ioff()  # this turns python interactive mode off to return to normal functions 
# pylab.plt.imshow(A,cmap='Greys',interpolation='nearest')
# pylab.plt.show()
# 
# ############################################################
#

#     setup cyclic  boundry conditions

# in fact in python you dont need to declare a third end of the array 
# if it is a np.array it will have a negative boundary and you only need to declare a positive boundry 
# for it to loop over itself

 bnd1 = np.array([np.arange(asize),np.arange(asize)])
 bnd2 = np.array([np.arange(bsize),np.arange(bsize)])
 bnd1 = bnd1.reshape(asize*2)
 bnd2 = bnd2.reshape(bsize*2)

# begining monte carlo sequence
  
# set parameters
 ngomax=asize*bsize
 iaccept1=0
 iaccept2=0
 ireject=0

# begin monte carlo loop 

 for mgo in range(ncycles):
  print 'mgo = %i %i %i %i' %(mgo,iaccept1,iaccept2,ireject)
  temp=1.0

# criteria for simulated annealling, reset counters 
  if(mgo>aneal): temp=0.1
  iaccept1=0
  iaccept2=0
  ireject=0
  
 #       dont need this bit 
 #       if(mgo.eq.1) call total(mgo)
 #       if(mgo.eq.1) call corr_abc(iocc,alpha,beta,side,ralpha,rbeta,rside)
  for ngo in range(ngomax):
  
  #  pick a block
   xx = np.random.random_integers(0,asize-1)
   yy = np.random.random_integers(0,bsize-1)

   e1 = energy(J,A,xx,yy,bnd1,bnd2)

  #  flip the value and calculate the energy 
   A[xx,yy] = -A[xx,yy]
  
   e2 = energy(J,A,xx,yy,bnd1,bnd2)

   dele=e2-e1
   rndvar = np.random.rand(1)
   
   if(dele < 0.0): 
     iaccept1 += 1
   elif(rndvar <= np.exp(-dele/temp)): 
     iaccept2 += 1
   else:
     ireject += 1
  #  swapback 
     A[xx,yy] = -A[xx,yy]

#  view each cycle
#  pylab.plt.imshow(A,cmap='Greys',interpolation='nearest')
#  pylab.plt.draw()
#  time.sleep(0.01)

 return A

def plot_grid(A):
# final result
 pylab.plt.ioff()  # this turns python interactive mode off to return to normal functions 
 pylab.plt.imshow(A,cmap='Greys',interpolation='nearest')
# pylab.plt.imshow(A,cmap='Accent',interpolation='nearest')
# pylab.plt.draw()  # just updates the plot and will turn off automatically 
# pylab.plt.show()  # just updates the plot and will turn off automatically 
 pylab.plt.savefig('occ_ising.png')  # just updates the plot and will turn off automatically 
# time.sleep(0.05)


def main():

 asize = 3
 bsize = 6
 L =[asize,bsize]

 Ja =  0.0
 Ja2 =  0.0
 Jb =  -0.1
 J = [Ja,Ja2,Jb]

 ncycles = 50
 aneal = 10
 temp = 0.1
 P = [ncycles,aneal,temp]

 A = ising(L,J,P)    # returns the occupacy square grid 

# the last two parameters are the size of the correlation function you want calculated 
 correlation(A,asize,bsize,4,4)

# write out the occ file for ZMC
# you need to dive into this routine each time and manually map the square grid to specific locations in the unit cell 
 cella = 3
 cellb = 3
 cellc = 3 
 occupancy(A,asize,bsize,cella,cellb,cellc)
 


if __name__ == "__main__":
    main()

############
