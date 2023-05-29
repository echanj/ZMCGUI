#!/usr/bin/perl 

use Math::Trig;
 
 my $pi = 3.1415926535897932384626433832795 ; 

############################
#
# shortened script extracted from make trial which only creates an editited set of diffuse input files 
# from the original 
#
# modify this script for differnt experiments by changing the content of the input files 
#
# ./edit_diffuse_input.pl trialno nframes sim1 sim2 sim3 lotsize1 lotsize2 lotsize3 nlots sub_ave_att 
#
# not that the simualtion size can be very different to the lotsize enableing larger simulation than lot 
#
##############################
print " to run type \n" ; 
print " ./edit_diffuse_input.pl trialno nframes sim1 sim2 sim3 lotsize1 lotsize2 lotsize3 nlots sub_ave_latt\n" ;  

#  my $rootname = $ARGV[0] ;
  my $trialno = $ARGV[0] ;
  my $nframes = $ARGV[1] ; # runtime in femtoseconds

  my $sim1 = $ARGV[2]; 
  my $sim2 = $ARGV[3]; 
  my $sim3 = $ARGV[4]; 

  my $lot1 = $ARGV[5]; 
  my $lot2 = $ARGV[6]; 
  my $lot3 = $ARGV[7]; 

  my $nlots = $ARGV[8];
  my $sub_ave_latt = $ARGV[9];

  my $frames = $nframes ;

########################################
# now edit the diffuse planes 
#################################################################
#
my @planes = ('h0l','0kl','hk0') ;
 
$simsize1 = $frames*$sim1;
#
foreach $plane (@planes) {

  open(diffin, "diffuse\_$plane.in") ;
  open(diffout, ">","diffuse\_$plane\_$trialno.in") ;

  while(<diffin>) { 
  if (/simulation_size/){ print diffout  "$simsize1 $sim2 $sim3  ! simulation_size\n" }
  elsif (/Lot size/){ print diffout "$lot1 $lot2 $lot3           ! lot_size\n" }
  elsif (/Number of lots/){ print diffout "$nlots                ! number of lots\n" }
  elsif (/Subtract average lattice/){ print diffout "$sub_ave_latt                ! Subtract average lattice\n" }
  else{print diffout $_ ;} 
  } ; 

  close (diffin) ;
  close (diffout) ;
    }

#
################################################################
#
#


