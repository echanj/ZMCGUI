#!/usr/bin/perl

# merges 2 mol2 files together they should be the same space group 
# make sure you have all the molecules in the order you want first 
# there is another program ive made for doing this  

#  $ program6_12 myfile1 myfile2
#  sets @ARGV to the list
#  ("myfile1", "myfile2")
#  $var = $ARGV[0];
#  assigns the first element of @ARGV to the scalar variable $var.



# @<TRIPOS>ATOM
#      1 C21      5.8059  -0.6904   2.9141   C.3       1 RES1   0.0000


$section1 = "@<TRIPOS>MOLECULE";
$section2 = "@<TRIPOS>ATOM";
$section3 = "@<TRIPOS>BOND";
$section4 = "@<TRIPOS>SUBSTRUCTURE";
$section5 = "@<TRIPOS>CRYSIN";
$count = 0 ;

open(indata1, $ARGV[0] );
open(indata2, $ARGV[1] );
open(head_1, ">head_1" );
open(tail_1, ">tail_1" );
open(head_2, ">head_2" );
open(tail_2, ">tail_2" );

while(<indata1>) {
#  next if (/$section1/); this will only single out that line
  if (/^\s*\d/) { 
         @h = split ;
         $num_atoms_1 = $h[0];
         $num_bonds_1 = $h[1];
         $num_subst_1 = $h[2];
                }
  print head_1 ; 
#  next unless (/$section2/);
  last if (/$section2/);
}

my (@num, @atom, @x, @y ,@z, @type, @mol, @res, @dum);
# @num = (); # initialise the arrays 
# @atom = ();
# @x = ();
# @y = ();
# @z = ();
# @type = ();
# @mol = ();
# @res = ();
# @dum = ();

# print " \n";
# close(indata1);

# open(indata1, $ARGV[0] );
while(<indata1>) {
#  next unless (/$section2/);
  if  (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @num,  $field[0];
  push @atom, $field[1];
  push @x ,   $field[2];
  push @y ,   $field[3];
  push @z ,   $field[4];
  push @type, $field[5];
  push @mol , $field[6];
  push @res , $field[7];
  push @dum , $field[8];
                 }
   last if (/$section3/) ; 
 }

# print " $num[3] \n";
# print "$num[0] $atom[0] $x[0] $y[0] $z[0] $type[0] $mol[0] $res[0] $dum[0]\n";
# sprintf "%.3f";
# sprintf "%3d";
# foreach my $i (@num) { 
#               print "$num[$i] $atom[$i] $x[$i] $y[$i] $z[$i] $type[$i] $mol[$i] $res[$i] $dum[$i]\n";
# }
# print "$#num \n" 

# for ($i=0; $i <= $#num ; $i++) {

#               print "$num[$i] $atom[$i] $x[$i] $y[$i] $z[$i] $type[$i] $mol[$i] $res[$i] $dum[$i]\n";
#}

my(@bond_id, @origin, @target, @bond_type);
 while(<indata1>) {
  if (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @bond_id,  $field[0];
  push @origin,  $field[1];
  push @target,  $field[2];
  push @bond_type ,  $field[3];
            }
   last if (/$section4/);
 }

# for ($i=0; $i <= $#bond_id ; $i++) {

#               print "$bond_id[$i] $origin[$i] $target[$i] $bond_type[$i] \n";
#}

my(@subst_id, @subst_name, @root_atom);
 while(<indata1>) {
  if (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @subst_id,  $field[0];
  push @subst_name,  $field[1];
  push @root_atom,  $field[2];
            }
   last if (/$section5/);
 }

#  for ($i=0; $i <= $#subst_id ; $i++) {
#                print "$subst_id[$i] $subst_name[$i] $root_atom[$i] \n";
# }


while(<indata1>) {
  print tail_1 ; 
}
close(indata1);


# now slurp data from second file 
 
 
 while(<indata2>) {
   if (/^\s*\d/) { 
          @h = split ;
          $num_atoms_2 = $h[0];
          $num_bonds_2 = $h[1];
          $num_subst_2 = $h[2];
                 }
   print head_2 ; 
   last if (/$section2/);
 }



my (@num2, @atom2, @x2, @y2 ,@z2, @type2, @mol2, @res2, @dum2);
while(<indata2>) {
  if  (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @num2,  $field[0];
  push @atom2, $field[1];
  push @x2 ,   $field[2];
  push @y2 ,   $field[3];
  push @z2 ,   $field[4];
  push @type2, $field[5];
  push @mol2 , $field[6];
  push @res2 , $field[7];
  push @dum2 , $field[8];
                 }
   last if (/$section3/) ; 
 }

# print " $num[3] \n";
# print "$num[0] $atom[0] $x[0] $y[0] $z[0] $type[0] $mol[0] $res[0] $dum[0]\n";
# sprintf "%.3f";
# sprintf "%3d";
# foreach my $i (@num) { 
#               print "$num[$i] $atom[$i] $x[$i] $y[$i] $z[$i] $type[$i] $mol[$i] $res[$i] $dum[$i]\n";
# }
# print "$#num \n" 

# for ($i=0; $i <= $#num2 ; $i++) {
#                print "$num2[$i] $atom2[$i] $x2[$i] $y2[$i] $z2[$i] $type2[$i] $mol2[$i] $res2[$i] $dum2[$i]\n";
# }

my(@bond_id2, @origin2, @target2, @bond_type2);
 while(<indata2>) {
  if (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @bond_id2,  $field[0];
  push @origin2,  $field[1];
  push @target2,  $field[2];
  push @bond_type2 ,  $field[3];
            }
   last if (/$section4/);
 }

# for ($i=0; $i <= $#bond_id ; $i++) {

# print "$bond_id[$i] $origin[$i] $target[$i] $bond_type[$i] \n";
#}

my(@subst_id2, @subst_name2, @root_atom2);
 while(<indata2>) {
  if (/^\s*\d/) {
#  print $_;
   chop; 
   s/^\s+//;
   @field = split; 
  push @subst_id2,  $field[0];
  push @subst_name2,  $field[1];
  push @root_atom2,  $field[2];
            }
   last if (/$section5/);
 }

#  for ($i=0; $i <= $#subst_id2 ; $i++) {
#               print "$subst_id2[$i] $subst_name2[$i] $root_atom2[$i] \n";
#  }


while(<indata2>) {
  print tail_2 ; 
}
close(indata2);



open(head_1, "head_1" );
open(tail_1, "tail_1" );
# 
$num_atoms = $num_atoms_1 + $num_atoms_2;
$num_bonds = $num_bonds_1 + $num_bonds_2;
$num_subst = $num_subst_1 + $num_subst_2;
# 
 while(<head_1>) {
   if (/^\s*\d/) { 
       print "   $num_atoms   $num_bonds   $num_subst \n" ;
  }else{ print ; } 
#   last if (/$section2/);
 }


for ($i=0; $i <= $#num ; $i++) {
                        if ($atom[$i] =~ /\d/){ 
                        $atom_label = "$atom[$i]" ;
                        }else{
                        $atom_label = "$atom[$i]$num[$i]" ;}
                print "   $num[$i] $atom_label $x[$i] $y[$i] $z[$i] $type[$i] $mol[$i] RES$mol[$i] $dum[$i]\n";
}

# print $#num, "\n";
# print $mol[$#num], "\n";

for ($i=0; $i <= $#num2 ; $i++) {
             $new_num2 = $num2[$i]+$num_atoms_1;
             $new_mol2 = $mol[$#num]+1;
             $new_res2 = "RES$new_mol2";
             $atom_label = "$atom2[$i]$num2[$i]" ;
                print "   $new_num2 $atom_label $x2[$i] $y2[$i] $z2[$i] $type2[$i] $new_mol2 $new_res2 $dum2[$i]\n";
}
#needto get rid of any whitespce after each section header 
chomp($section3);
print "$section3\n" ;

for ($i=0; $i <= $#bond_id ; $i++) {
print "  $bond_id[$i] $origin[$i] $target[$i] $bond_type[$i] \n";
}

for ($i=0; $i <= $#bond_id2 ; $i++) {
     $new_bond_id2 = $bond_id2[$i]+$num_bonds_1;
     $new_origin2 = $origin2[$i]+$num_atoms_1;
     $new_target2 = $target2[$i]+$num_atoms_1;
print "  $new_bond_id2 $new_origin2 $new_target2 $bond_type2[$i] \n";
}

chomp($section4);
print "$section4\n" ;

for ($i=0; $i <= $#subst_id ; $i++) {
             print "   $subst_id[$i] $subst_name[$i] $root_atom[$i] GROUP  0 ****  ****  0 \n";
}

for ($i=0; $i <= $#subst_id2 ; $i++) {
              $new_subst_id2 = $subst_id2[$i]+ $subst_id[$#subst_id]   ;
              $new_subst_name2 = "RES$new_subst_id2"   ;
              $new_root_atom2 =  $root_atom2[$i]+ $num_atoms_1   ;
             print "   $new_subst_id2 $new_subst_name2 $new_root_atom2 GROUP  0 ****  ****  0 \n";
}

chomp($section5);
print "$section5\n" ;
while(<tail_1>){print;}

unlink(head_1);
unlink(tail_1);
# 
unlink(head_2);
unlink(tail_2);
# 
# 
# 
