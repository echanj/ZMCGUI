#!/bin/bash

# simple bash script to run calculate diffraction images  
# you should edit this file prior 
# then use it to create the template by running diffuse input gen

#
# remove files for previous run to prevent conflicts 
# 

rm input_diffuse*
rm INPUT_*
rm diffuse_*.in

# currently setup for lammps2diff 

cat > INPUT_H0L << EOF
aspirin_formI_h0l                      ! title 
diffuse_h0l.in                                 ! file_name
11.416   6.598  11.483  90  95.600  90 ! unit cell  
0.7000 25.5                            ! lambda  theta_max 
0.0 0.0 0.0                            ! origin
1.0 0.0 0.0                            ! vertical axis
0.0 0.0 1.0                            ! horizontal axis
EOF

cat > INPUT_HK0 << EOF2
aspirin_formI_hk0                      ! title 
diffuse_hk0.in                                 ! file_name
11.416   6.598  11.483  90  95.600  90 ! unit cell  
0.7000 25.5                            ! lambda  theta_max 
0.0 0.0 0.0                            ! origin
0.0 1.0 0.0                            ! vertical axis
1.0 0.0 0.0                            ! horizontal axis
EOF2

cat > INPUT_0KL << EOF3
aspirin_formI_0kl                      ! title 
diffuse_0kl.in                                 ! file_name
11.416   6.598  11.483  90  95.600  90 ! unit cell  
0.7000 25.5                            ! lambda  theta_max 
0.0 0.0 0.0                            ! origin
0.0 0.0 1.0                            ! vertical axis
0.0 1.0 0.0                            ! horizontal axis
EOF3

# $cell_lengths 90 68.491 90  ! unit cell  

# awk '{ sub(/8                              ! Number of atom sites per cell/, "104"); print }' diffuse_test.in > diffuse_h0l.in

cat > input_diffuse_h0l << EOF
diffuse_h0l.in
h0l.bin
EOF

cat > input_diffuse_hk0 << EOF
diffuse_hk0.in
hk0.bin
EOF

cat > input_diffuse_0kl << EOF
diffuse_0kl.in
0kl.bin
EOF

# ../ZMC --diffuse MC_params_mol.inp aspirin

diffuse_input_generator < INPUT_H0L
diffuse_input_generator < INPUT_HK0
diffuse_input_generator < INPUT_0KL

echo 'this is an instance of find_lot_size calculation'  
echo 'it can be used as a guide to help figure out what lot size parameters to use'
echo 'it is currently tuned to a crystal size of 48 48 48 unitcells'
find_lot_sizes.pl  11.416   6.598  11.483   48  48  48  15 

cat > run_DZMC.sh << EOFF
#!/bin/bash

# make sure you correctly edit input files 

ZMC --crystal --diffuse form1_ZMC.inp form1


DZMC form1.diffuse < input_diffuse_h0l  
~/TRW_bundle/bin/bin2gray  --twofold --norm=20000 h0l.bin 
mv h0l.pgm pgm/h0l.pgm
rm h0l.bin 

DZMC form1.diffuse < input_diffuse_hk0  
~/TRW_bundle/bin/bin2gray --hmirror -vmirror  --norm=20000 hk0.bin 
mv hk0.pgm pgm/hk0.pgm
rm hk0.bin 

DZMC form1.diffuse < input_diffuse_0kl  
~/TRW_bundle/bin/bin2gray --hmirror -vmirror  --norm=20000 0kl.bin 
mv 0kl.pgm pgm/0kl.pgm
rm 0kl.bin 

EOFF

chmod +x run_DZMC.sh 


##########################################################
# now a short script to change values which are identical in all the diffuse files
# on the sed line you can change the y to an n or an e 
# descriptions :
# line 4 :   !simulation size 
# line 11 :  !lot size 
# line 12 :  !number of lots 
# line 13 :  !number of atoms per cell 
# line 14 :  !number of atom types  
# line 15 :  !subtract average lattice 
###########################################

file="diffuse_*.in"

for i in $file  
do 
temp=$i.temp
temp2=$i.temp2
temp3=$i.temp3
temp4=$i.temp4
temp5=$i.temp5
temp6=$i.temp6
sed '4 s/48 48 48/20 20 20/' $i > $temp
sed '11 s/8,15,8/10,10,10/' $temp > $temp2
sed '12 s/1/1/' $temp2 > $temp3
sed '13 s/8/84/' $temp3 > $temp4
sed '14 s/2/2/' $temp4 > $temp5
sed '15 s/n/n/' $temp5 > $temp6
cp $temp6 $i
rm $temp
rm $temp2
rm $temp3
rm $temp4
rm $temp5
rm $temp6
done

# delete unwanted files to save space 
# rm *.bin 
# rm input_diffuse*


#  
#



