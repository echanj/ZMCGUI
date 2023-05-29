#!/usr/bin/perl

# enter the cell parameters and crystal dimensions and 
# program will produce a table for a distribuion of possible lot size values that 
# should be of equal dimensions  
# 
# routine for finding the maximum value 
sub max { @_ = sort { $a <=> $b } @_; pop }

# number of iterations
$num =  $ARGV[6] ;

$crya =  $ARGV[3]  ;
$cryb =  $ARGV[4] ;
$cryc =  $ARGV[5] ; 

$cella = $ARGV[0] ;
$cellb = $ARGV[1] ;
$cellc = $ARGV[2] ;

# feed cell constants in to an array 
@cell = ($cella,$cellb,$cellc);

$cell_max = max(@cell);

# print @sorted_cell;

# $cell_min = $sorted_cell[0] 
# $cell_mid = $sorted_cell[1] 
# $cell_max = $sorted_cell[2] 

print "$cell_max \n"; 

$c1= $cell_max/$cella; 
$c2= $cell_max/$cellb; 
$c3= $cell_max/$cellc; 

$count=1;

$crysize=$crya*$cryb*$cryc;

  for ($i=1; $i <= $num ; $i++) {

   $lota = $c1*$count ; 
   $lotb = $c2*$count ;
   $lotc = $c3*$count ;

  $nlot = $crysize/($lota*$lotb*$lotc) ;

print "lot_abc = $lota $lotb $lotc num_lots = $nlot \n";
$count++;

}


#
