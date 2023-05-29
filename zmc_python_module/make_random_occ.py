#!/opt/local/bin/python

# run by typing 
#
# >>python make_random_occ.py cellsize loc zocc mocc
#
# where 
# 
# cellsize = the size of the crystal in this case it is cubic  
# loc = number of locations 
# zocc = the number of zmatrices
# mocc = the numbe rof molecules per location
#

import numpy as np

from sys import argv 

script, cellsize, loc, zocc, mocc = argv


print "run by typing" 
print "python make_random_occ.py cellsize loc zocc mocc"


# convert the input to and integer 
cella = int(cellsize) 
cellb = int(cellsize)
cellc = int(cellsize) 
loc = int(loc) 
zocc = int(zocc) 
mocc = int(mocc) 

# use np.random.random_integers() in order to get a number from 1-zocc

for a in range(cella):
 for b in range(cellb):
  for c in range(cellc):
   for l in range(loc):
    z = np.random.random_integers(1,zocc,1)
    m = np.random.random_integers(1,mocc,1)
    print a+1, b+1, c+1, l+1, z[0], m[0]




