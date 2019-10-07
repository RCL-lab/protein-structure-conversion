####################################################################################################
## FileName:    Extract-Sidechain-internalCordinate.tcl
## Description: This module extract the internal coordinate of both protein Sidechain
##				Load protein.pdb input in VMD 
##				source Extract-Sidechain-internalCordinate.tcl
##				type the number of residue of your protein as input
##				The output are 5 files "BOND_SIDE.txt", "ANGLE_SIDE.txt", "imprANGLE_SIDE.txt", "DIHD_SIDE.txt", "resname.txt"
## 				
## Author: Mahsa Bayati
##
## MIT License
## Copyright (c) 2019 RCL-lab
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
####################################################################################################

puts "Enter the protein residue No.: "
gets stdin NoRes




set selatomN [[atomselect top "name N"] get index ]
set selatomCA [[atomselect top "name CA"] get index ]
set selatomC [[atomselect top "name C"] get index ]

set ALA "ALA"
set ARG "ARG"
set ASN "ASN"
set ASP "ASP"
set CYS "CYS"
set GLU "GLU"
set GLN "GLN"
set GLY "GLY"
set HSD "HSD"
set ILE "ILE"
set LEU "LEU"
set LYS "LYS"
set MET "MET"
set PHE "PHE"
set PRO "PRO"
set SER "SER"
set THR "THR"
set TRP "TRP"
set TYR "TYR"
set VAL "VAL"
set HIS "HIS"


set bondmeasure [open "BOND_SIDE.txt" "w"]
set anglemeasure [open "ANGLE_SIDE.txt" "w"]
set impranglemeasure [open "imprANGLE_SIDE.txt" "w"]
set dihdmeasure [open "DIHD_SIDE.txt" "w"]
set resi [open "resname.txt" "w"]
for {set i 2} {$i < $NoRes} {incr i} {
	set side_i [atomselect top "resid $i"];
	set s_TF [$side_i get { sidechain}];
	set s_name  [$side_i get { name}];
	set s_index [$side_i get {index}];
	set res_name [$side_i get {resname}];

	if {[string compare [lindex $res_name 0] ALA]==0} {
	puts   $resi "$ALA";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set HB3 [lindex $s_index $c_atom];

	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]
	
    set bCACB [measure bond [list $CA $CB]];  
	puts   $bondmeasure "$bCACB";
	
	set bCBHB1 [measure bond [list $CB $HB1]];  
	puts   $bondmeasure "$bCBHB1";
	
    set bCBHB2 [measure bond [list $CB $HB2]];  
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCBHB3 [measure bond [list $CB $HB3]];  
	puts   $bondmeasure "$bCBHB3";
	
	for {set j 0} {$j <19} {incr j} {
	puts   $bondmeasure "0";
    }
	
	set aNCACB [measure angle [list $N $CA $CB]];   
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBHB3 [measure angle [list $CA $CB $HB3]];  
	puts   $anglemeasure "$aCACBHB3";
	
	for {set j 0} {$j <19} {incr j} {
	puts   $anglemeasure "0";
    }
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; ; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dNCACBHB3 [measure dihed [list $N $CA $CB $HB3]]; 
	puts   $dihdmeasure "$dNCACBHB3";
	
	for {set j 0} {$j <19} {incr j} {
	puts   $dihdmeasure "0";
    }
    }

	if {[string compare [lindex $res_name 0] THR]==0} {
	puts   $resi "$THR";
	set c_atom 4
	set CB [lindex $s_index $c_atom];

	incr c_atom;
	set HB [lindex $s_index $c_atom];
	incr c_atom;
	set OG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
   
	
	incr c_atom;
	set CG2 [lindex $s_index $c_atom];
   	incr c_atom;
	set HG21 [lindex $s_index $c_atom];
	incr c_atom;
	set HG22 [lindex $s_index $c_atom];
	incr c_atom;
	set HG23 [lindex $s_index $c_atom];

	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB [measure bond [list $CB $HB]];  
	puts   $bondmeasure "$bCBHB";
	
    set bCBOG1 [measure bond [list $CB $OG1]];  
	puts   $bondmeasure "$bCBOG1";
	
	set bOG1HG1 [measure bond [list $OG1 $HG1]];  
	puts   $bondmeasure "$bOG1HG1";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCBCG2 [measure bond [list $CB $CG2]];  
	puts   $bondmeasure "$bCBCG2";
	
	set bCG2HG21 [measure bond [list $CG2 $HG21]];  
	puts   $bondmeasure "$bCG2HG21";
	
	set bCG2HG22 [measure bond [list $CG2 $HG22]];   
	puts   $bondmeasure "$bCG2HG22";

	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCG2HG23 [measure bond [list $CG2 $HG23]];  
	puts   $bondmeasure "$bCG2HG23";
	
	for {set j 0} {$j <13} {incr j} {
	puts   $bondmeasure "0";
	}
	set aNCACB [measure angle [list $N $CA $CB]];   
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB [measure angle [list $CA $CB $HB]];  
	puts   $anglemeasure "$aCACBHB";	
	
	set aCACBOG1 [measure angle [list $CA $CB $OG1]];  
	puts   $anglemeasure "$aCACBOG1";
	
	set aCBOG1HG1 [measure angle [list $CB $OG1 $HG1]];  
	puts   $anglemeasure "$aCBOG1HG1";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG2 [measure angle [list $CA $CB $CG2]];  
	puts   $anglemeasure "$aCACBCG2";
	
	set aCBCG2HG21 [measure angle [list $CB $CG2 $HG21]];   
	puts   $anglemeasure "$aCBCG2HG21";
	
	set aCBCG2HG22 [measure angle [list $CB $CG2 $HG22]];  
	puts   $anglemeasure "$aCBCG2HG22";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCG2HG23 [measure angle [list $CB $CG2 $HG23]];  
	puts   $anglemeasure "$aCBCG2HG23";
   
    for {set j 0} {$j <13} {incr j} {
	puts   $anglemeasure "0";
    }
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB [measure dihed [list $N $CA $CB $HB]]; 
	puts   $dihdmeasure "$dNCACBHB";
	
	
	set dNCACBOG1 [measure dihed [list $N $CA $CB $OG1]]; 
	puts   $dihdmeasure "$dNCACBOG1";
	
	set dCACBOG1HG1 [measure dihed [list $CA $CB $OG1 $HG1]]; 
	puts   $dihdmeasure "$dCACBOG1HG1";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dNCACBCG2 [measure dihed [list $N $CA $CB $CG2]]; 
	puts   $dihdmeasure "$dNCACBCG2";
		
	set dCACBCG2HG21 [measure dihed [list $CA $CB $CG2 $HG21]];  
	puts   $dihdmeasure "$dCACBCG2HG21";
	
	set dCACBCG2HG22 [measure dihed [list $CA $CB $CG2 $HG22]]; 
	puts   $dihdmeasure "$dCACBCG2HG22";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCG2HG23 [measure dihed [list $CA $CB $CG2 $HG23]]; 
	puts   $dihdmeasure "$dCACBCG2HG23";
	for {set j 0} {$j <13} {incr j} {
	puts   $dihdmeasure "0";
    }

    }

	if {[string compare [lindex $res_name 0] LEU]==0} {
	puts   $resi "$LEU";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set HG [lindex $s_index $c_atom];
	incr c_atom;
	set CD1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD11 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD12 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD13 [lindex $s_index $c_atom];
   	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
   	incr c_atom;	
	set HD21 [lindex $s_index $c_atom];
	incr c_atom;
	set HD22 [lindex $s_index $c_atom];
	incr c_atom;
	set HD23 [lindex $s_index $c_atom];
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG [measure bond [list $CG $HG]];  
	puts   $bondmeasure "$bCGHG";
	
	set bCGCD1 [measure bond [list $CG $CD1]];  
	puts   $bondmeasure "$bCGCD1";
	
	set bCD1HD11 [measure bond [list $CD1 $HD11]];  
	puts   $bondmeasure "$bCD1HD11";
	
	set bCD1HD12 [measure bond [list $CD1 $HD12]];  
	puts   $bondmeasure "$bCD1HD12";
	
	set bCD1HD13 [measure bond [list $CD1 $HD13]];  
	puts   $bondmeasure "$bCD1HD13";
	
	set bCGCD2 [measure bond [list $CG $CD2]];  
	puts   $bondmeasure "$bCGCD2";

	set bCD2HD21 [measure bond [list $CD2 $HD21]];  
	puts   $bondmeasure "$bCD2HD21";
	
	set bCD2HD22 [measure bond [list $CD2 $HD22]];  
	puts   $bondmeasure "$bCD2HD22";
	
	set bCD2HD23 [measure bond [list $CD2 $HD23]];   
	puts   $bondmeasure "$bCD2HD23";
	
	for {set j 0} {$j <10} {incr j} {
	puts   $bondmeasure "0";
    }
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG [measure angle [list $CB $CG $HG]];  
	puts   $anglemeasure "$aCBCGHG";
	
	set aCBCGCD1 [measure angle [list $CB $CG $CD1]];  
	puts   $anglemeasure "$aCBCGCD1";
	
	set aCGCD1HD11 [measure angle [list $CG $CD1 $HD11]];  
	puts   $anglemeasure "$aCGCD1HD11";
	
	set aCGCD1HD12 [measure angle [list $CG $CD1 $HD12]];  
	puts   $anglemeasure "$aCGCD1HD12";
	
	set aCGCD1HD13 [measure angle [list $CG $CD1 $HD13]];  
	puts   $anglemeasure "$aCGCD1HD13";
	
	set aCBCGCD2 [measure angle [list $CB $CG $CD2]];  
	puts   $anglemeasure "$aCBCGCD2";
	
	set aCGCD2HD21 [measure angle [list $CG $CD2 $HD21]];  
	puts   $anglemeasure "$aCGCD2HD21";
	
	set aCGCD2HD22 [measure angle [list $CG $CD2 $HD22]];  
	puts   $anglemeasure "$aCGCD2HD22";
	
	set aCGCD2HD23 [measure angle [list $CG $CD2 $HD23]];  
	puts   $anglemeasure "$aCGCD2HD23";
	
	for {set j 0} {$j <10} {incr j} {
    puts   $anglemeasure "0";
    }
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	
	set dCACBCGHG [measure dihed [list $CA $CB $CG $HG]]; 
	puts   $dihdmeasure "$dCACBCGHG";
		
	set dCACBCGCD1 [measure dihed [list $CA $CB $CG $CD1]]; 
	puts   $dihdmeasure "$dCACBCGCD1";
	
	set dCBCGCD1HD11 [measure dihed [list $CB $CG $CD1 $HD11]]; 
	puts   $dihdmeasure "$dCBCGCD1HD11";
	
	set dCBCGCD1HD12 [measure dihed [list $CB $CG $CD1 $HD12]];  
	puts   $dihdmeasure "$dCBCGCD1HD12";
	
	set dCBCGCD1HD13 [measure dihed [list $CB $CG $CD1 $HD13]]; 
	puts   $dihdmeasure "$dCBCGCD1HD13";
	
	set dCACBCGCD2 [measure dihed [list $CA $CB $CG $CD2]]; 
	puts   $dihdmeasure "$dCACBCGCD2";
	
	set dCBCGCD2HD21 [measure dihed [list $CB $CG $CD2 $HD21]];  
	puts   $dihdmeasure "$dCBCGCD2HD21";
	
	set dCBCGCD2HD22 [measure dihed [list $CB $CG $CD2 $HD22]];  
	puts   $dihdmeasure "$dCBCGCD2HD22";
	
	set dCBCGCD2HD23 [measure dihed [list $CB $CG $CD2 $HD23]];  
	puts   $dihdmeasure "$dCBCGCD2HD23";
	
	for {set j 0} {$j <10} {incr j} {
	puts   $dihdmeasure "0";
    }

    }	

	if {[string compare [lindex $res_name 0] LYS]==0} {
	puts   $resi "$LYS";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];
	incr c_atom;
	set CE [lindex $s_index $c_atom];
   	incr c_atom;
	set NZ [lindex $s_index $c_atom];
	incr c_atom;
	incr c_atom;
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG2 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
   	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE2 [lindex $s_index $c_atom];
	incr c_atom;
	set HZ1 [lindex $s_index $c_atom];
	incr c_atom;
	set HZ2 [lindex $s_index $c_atom];
	incr c_atom;
	set HZ3 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG1 [measure bond [list $CG $HG1]];  
	puts   $bondmeasure "$bCGHG1";
	
	set bCGHG2 [measure bond [list $CG $HG2]];  
	puts   $bondmeasure "$bCGHG2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGCD [measure bond [list $CG $CD]];  
	puts   $bondmeasure "$bCGCD";
	
	set bCDHD1 [measure bond [list $CD $HD1]];  
	puts   $bondmeasure "$bCDHD1";
	set bCDHD2 [measure bond [list $CD $HD2]];  
	puts   $bondmeasure "$bCDHD2";
	
	set bCDCE [measure bond [list $CD $CE]];  
	puts   $bondmeasure "$bCDCE";
	
	set bCEHE1 [measure bond [list $CE $HE1]];  
	puts   $bondmeasure "$bCEHE1";
	set bCEHE2 [measure bond [list $CE $HE2]];  
	puts   $bondmeasure "$bCEHE2";
	
	
	set bCENZ [measure bond [list $CE $NZ]];  
	puts   $bondmeasure "$bCENZ";
	
	set bNZHZ1 [measure bond [list $NZ $HZ1]];  
	puts   $bondmeasure "$bNZHZ1";	
	set bNZHZ2 [measure bond [list $NZ $HZ2]];  
	puts   $bondmeasure "$bNZHZ2";
		
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bNZHZ3 [measure bond [list $NZ $HZ3]];  
	puts   $bondmeasure "$bNZHZ3";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG1 [measure angle [list $CB $CG $HG1]];  
	puts   $anglemeasure "$aCBCGHG1";
	
	set aCBCGHG2 [measure angle [list $CB $CG $HG2]];  
	puts   $anglemeasure "$aCBCGHG2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGCD [measure angle [list $CB $CG $CD]];  
	puts   $anglemeasure "$aCBCGCD";
		
	set aCGCDHD1 [measure angle [list $CG $CD $HD1]];  
	puts   $anglemeasure "$aCGCDHD1";
	
	set aCGCDHD2 [measure angle [list $CG $CD $HD2]];  
	puts   $anglemeasure "$aCGCDHD2";
	
	set aCGCDCE [measure angle [list $CG $CD $CE]];  
	puts   $anglemeasure "$aCGCDCE";
	
	set aCDCEHE1 [measure angle [list $CD $CE $HE1]];  
	puts   $anglemeasure "$aCDCEHE1";
	
	set aCDCEHE2 [measure angle [list $CD $CE $HE2]];  
	puts   $anglemeasure "$aCDCEHE2";
	
	set aCDCENZ [measure angle [list $CD $CE $NZ]];  
	puts   $anglemeasure "$aCDCENZ";
	
	set aCENZHZ1 [measure angle [list $CE $NZ $HZ1]];  
	puts   $anglemeasure "$aCENZHZ1";
	set aCENZHZ2 [measure angle [list $CE $NZ $HZ2]];  
	puts   $anglemeasure "$aCENZHZ2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	set aCENZHZ3 [measure angle [list $CE $NZ $HZ3]];  
	puts   $anglemeasure "$aCENZHZ3";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";

	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGHG1 [measure dihed [list $CA $CB $CG $HG1]]; 
	puts   $dihdmeasure "$dCACBCGHG1";
	set dCACBCGHG2 [measure dihed [list $CA $CB $CG $HG2]]; 
	puts   $dihdmeasure "$dCACBCGHG2";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCGCD [measure dihed [list $CA $CB $CG $CD]]; 
	puts   $dihdmeasure "$dCACBCGCD";
	
	set dCBCGCDHD1 [measure dihed [list $CB $CG $CD $HD1]]; 
	puts   $dihdmeasure "$dCBCGCDHD1";
	set dCBCGCDHD2 [measure dihed [list $CB $CG $CD  $HD2]]; 
	puts   $dihdmeasure "$dCBCGCDHD2";
	
	set dCBCGCDCE [measure dihed [list $CB $CG $CD $CE]]; 
	puts   $dihdmeasure "$dCBCGCDCE";
	
	set dCGCDCEHE1 [measure dihed [list $CG $CD $CE $HE1]]; 
	puts   $dihdmeasure "$dCGCDCEHE1";
	set dCGCDCEHE2 [measure dihed [list $CG $CD $CE $HE2]]; 
	puts   $dihdmeasure "$dCGCDCEHE2";
	
	set dCGCDCENZ [measure dihed [list $CG $CD $CE $NZ]]; 
	puts   $dihdmeasure "$dCGCDCENZ";
	
	set dCDCENZHZ1 [measure dihed [list $CD $CE $NZ $HZ1]]; 
	puts   $dihdmeasure "$dCDCENZHZ1";	
	set dCDCENZHZ2 [measure dihed [list $CD $CE $NZ $HZ2]]; 
	puts   $dihdmeasure "$dCDCENZHZ2";	
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCDCENZHZ3 [measure dihed [list $CD $CE $NZ $HZ3]]; 
	puts   $dihdmeasure "$dCDCENZHZ3";	

	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
    }
		
	if {[string compare [lindex $res_name 0] GLU]==0} {
	puts   $resi "$GLU";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
    set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];
	incr c_atom;
	set OE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set OE2 [lindex $s_index $c_atom];
   	incr c_atom;
   	incr c_atom;
   	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG2 [lindex $s_index $c_atom];
	incr c_atom;

	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG1 [measure bond [list $CG $HG1]];  
	puts   $bondmeasure "$bCGHG1";
	
	set bCGHG2 [measure bond [list $CG $HG2]];  
	puts   $bondmeasure "$bCGHG2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGCD [measure bond [list $CG $CD]];  
	puts   $bondmeasure "$bCGCD";
	
	set bCDOE1 [measure bond [list $CD $OE1]];  
	puts   $bondmeasure "$bCDOE1";
	set bCDOE2 [measure bond [list $CD $OE2]];  
	puts   $bondmeasure "$bCDOE2";
	
	for {set j 0} {$j <11} {incr j} {
	puts   $bondmeasure "0";
    }
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG1 [measure angle [list $CB $CG $HG1]];  
	puts   $anglemeasure "$aCBCGHG1";
	
	set aCBCGHG2 [measure angle [list $CB $CG $HG2]];  
	puts   $anglemeasure "$aCBCGHG2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGCD [measure angle [list $CB $CG $CD]];  
	puts   $anglemeasure "$aCBCGCD";
		
	set aCGCDOE1 [measure angle [list $CG $CD $OE1]];  
	puts   $anglemeasure "$aCGCDOE1";
	
	set aCGCDOE2 [measure angle [list $CG $CD $OE2]];  
	puts   $anglemeasure "$aCGCDOE2";
	
	for {set j 0} {$j <11} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGHG1 [measure dihed [list $CA $CB $CG $HG1]]; 
	puts   $dihdmeasure "$dCACBCGHG1";
	set dCACBCGHG2 [measure dihed [list $CA $CB $CG $HG2]]; 
	puts   $dihdmeasure "$dCACBCGHG2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCGCD [measure dihed [list $CA $CB $CG $CD]]; 
	puts   $dihdmeasure "$dCACBCGCD";
	
	set dCBCGCDOE1 [measure dihed [list $CB $CG $CD  $OE1]]; 
	puts   $dihdmeasure "$dCBCGCDOE1";
	set dCBCGCDOE2 [measure dihed [list $CB $CG $CD $OE2]]; 
	puts   $dihdmeasure "$dCBCGCDOE2";

	for {set j 0} {$j <11} {incr j} {
	puts   $dihdmeasure "0";
    }

    }

	if {[string compare [lindex $res_name 0] GLY]==0} {
	puts   $resi "$GLY";
	set c_atom 4
	set HA2 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCAHA2 [measure bond [ list $CA $HA2 ]];
	puts  $bondmeasure "$bCAHA2"

    for {set j 0} {$j <25} {incr j} {
	puts   $bondmeasure "0";
    }
	
	set aNCAHA2 [measure angle [list $N $CA $HA2]];  
	puts   $anglemeasure "$aNCAHA2";
	
	for {set j 0} {$j <25} {incr j} {
	puts   $anglemeasure "0";
    }
	
	set aHA2CAC [measure angle [list $HA2 $CA $C]];  
	puts   $impranglemeasure "$aHA2CAC";
	
	set dCANHA2C [measure dihed [list $CA $N $HA2 $C]]; 
	puts   $dihdmeasure "$dCANHA2C";
	
	for {set j 0} {$j <25} {incr j} {
	puts   $dihdmeasure "0";
	}

    }

	if {[string compare [lindex $res_name 0] ASP]==0} {
	puts   $resi "$ASP";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set OD1 [lindex $s_index $c_atom];
	incr c_atom;
	set OD2 [lindex $s_index $c_atom];
	incr c_atom;
	incr c_atom;
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];



	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGOD1 [measure bond [list $CG $OD1]];  
	puts   $bondmeasure "$bCGOD1";
	
	set bCGOD2 [measure bond [list $CG $OD2]];  
	puts   $bondmeasure "$bCGOD2";
	
	for {set j 0} {$j <17} {incr j} {
	puts   $bondmeasure "0";
    }
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGOD1 [measure angle [list $CB $CG $OD1]];  
	puts   $anglemeasure "$aCBCGOD1";
	
	set aCBCGOD2 [measure angle [list $CB $CG $OD2]];  
	puts   $anglemeasure "$aCBCGOD2";
	
	for {set j 0} {$j <17} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGOD1 [measure dihed [list $CA $CB $CG $OD1]]; 
	puts   $dihdmeasure "$dCACBCGOD1";
	set dCACBCGOD2 [measure dihed [list $CA $CB $CG $OD2]]; 
	puts   $dihdmeasure "$dCACBCGOD2";
	
	for {set j 0} {$j <17} {incr j} {
	puts   $dihdmeasure "0";
    }

    }
	
	if {[string compare [lindex $res_name 0] CYS]==0} {
	puts   $resi "$CYS";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set SG [lindex $s_index $c_atom];
	incr c_atom;
	incr c_atom;
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];


	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBSG [measure bond [list $CB $SG]];  
	puts   $bondmeasure "$bCBSG";
	
	set bSGHG1 [measure bond [list $SG $HG1]];  
	puts   $bondmeasure "$bSGHG1";
	
for {set j 0} {$j <18} {incr j} {	
	puts   $bondmeasure "0";	
}
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBSG [measure angle [list $CA $CB $SG]];   
	puts   $anglemeasure "$aCACBSG";
	
	set aCBSGHG1 [measure angle [list $CB $SG $HG1]];  
	puts   $anglemeasure "$aCBSGHG1";
	
    for {set j 0} {$j <18} {incr j} {	
	puts   $anglemeasure "0";
    }
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBSG [measure dihed [list $N $CA $CB $SG]]; 
	puts   $dihdmeasure "$dNCACBSG";
		
	set dCACBSGHG1 [measure dihed [list $CA $CB $SG $HG1]]; 
	puts   $dihdmeasure "$dCACBSGHG1";

	for {set j 0} {$j <18} {incr j} {	
	puts   $dihdmeasure "0";
    }
    }
	if {[string compare [lindex $res_name 0] SER]==0} {
	puts   $resi "$SER";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set OG [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];


	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBOG [measure bond [list $CB $OG]];  
	puts   $bondmeasure "$bCBOG";
	
	set bOGHG1 [measure bond [list $OG $HG1]];  
	puts   $bondmeasure "$bOGHG1";
 for {set j 0} {$j <18} {incr j} {
	puts   $bondmeasure "0";	
 }
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBOG [measure angle [list $CA $CB $OG]];   
	puts   $anglemeasure "$aCACBOG";
	
	set aCBOGHG1 [measure angle [list $CB $OG $HG1]];  
	puts   $anglemeasure "$aCBOGHG1";
	
    for {set j 0} {$j <18} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBOG [measure dihed [list $N $CA $CB $OG]]; 
	puts   $dihdmeasure "$dNCACBOG";
		
	set dCACBOGHG1 [measure dihed [list $CA $CB $OG $HG1]]; 
	puts   $dihdmeasure "$dCACBOGHG1";

	for {set j 0} {$j <18} {incr j} {
	puts   $dihdmeasure "0";
	}
    }	
	if {[string compare [lindex $res_name 0] ILE]==0} {
	puts   $resi "$ILE";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set CG1 [lindex $s_index $c_atom];
	incr c_atom;
	set CG2 [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];
	incr c_atom;	
	incr c_atom;
	incr c_atom;
	set HB [lindex $s_index $c_atom];
	incr c_atom;
	set HG11 [lindex $s_index $c_atom];
	incr c_atom;
	set HG12 [lindex $s_index $c_atom];
	incr c_atom;
	set HG21 [lindex $s_index $c_atom];
	incr c_atom;
	set HG22 [lindex $s_index $c_atom];
	incr c_atom;
	set HG23 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
	incr c_atom;
	set HD3 [lindex $s_index $c_atom];

	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB [measure bond [list $CB $HB]];   
	puts   $bondmeasure "$bCBHB";
	
	set bCBCG2 [measure bond [list $CB $CG2]];   
	puts   $bondmeasure "$bCBCG2";
	
	set bCG2HG21 [measure bond [list $CG2 $HG21]];   
	puts   $bondmeasure "$bCG2HG21";
	
	set bCG2HG22 [measure bond [list $CG2 $HG22]];   
	puts   $bondmeasure "$bCG2HG22";
	
	set bCG2HG23 [measure bond [list $CG2 $HG23]];   
	puts   $bondmeasure "$bCG2HG23";
	
	set bCBCG1 [measure bond [list $CB $CG1]];   
	puts   $bondmeasure "$bCBCG1";
	
	set bCG1HG11 [measure bond [list $CG1 $HG11]];   
	puts   $bondmeasure "$bCG1HG11";
	
	set bCG1HG12 [measure bond [list $CG1 $HG12]];   
	puts   $bondmeasure "$bCG1HG12";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCG1CD [measure bond [list $CG1 $CD]];  
	puts   $bondmeasure "$bCG1CD";
	
	set bCDHD1 [measure bond [list $CD $HD1]];  
	puts   $bondmeasure "$bCDHD1";
	set bCDHD2 [measure bond [list $CD $HD2]];  
	puts   $bondmeasure "$bCDHD2";
	set bCDHD3 [measure bond [list $CD $HD3]];  
	puts   $bondmeasure "$bCDHD3";

	for {set j 0} {$j <10} {incr j} {
	puts   $bondmeasure "0";	
	}
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB [measure angle [list $CA $CB $HB]];  
	puts   $anglemeasure "$aCACBHB";	
	
	set aCACBCG2 [measure angle [list $CA $CB $CG2]];  
	puts   $anglemeasure "$aCACBCG2";
	set aCBCG2HG21 [measure angle [list $CB $CG2 $HG21]];  
	puts   $anglemeasure "$aCBCG2HG21";
	set aCBCG2HG22 [measure angle [list $CB $CG2 $HG22]];  
	puts   $anglemeasure "$aCBCG2HG22";
	set aCBCG2HG23 [measure angle [list $CB $CG2 $HG23]];  
	puts   $anglemeasure "$aCBCG2HG23";
	
	set aCACBCG1 [measure angle [list $CA $CB $CG1]];  
	puts   $anglemeasure "$aCACBCG1";
	set aCBCG1HG11 [measure angle [list $CB $CG1 $HG11]];  
	puts   $anglemeasure "$aCBCG1HG11";
	set aCBCG1HG12 [measure angle [list $CB $CG1 $HG12]];  
	puts   $anglemeasure "$aCBCG1HG12";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCG1CD [measure angle [list $CB $CG1 $CD]];  
	puts   $anglemeasure "$aCBCG1CD";
	
	set aCG1CDHD1 [measure angle [list $CG1 $CD $HD1]];  
	puts   $anglemeasure "$aCG1CDHD1";
	set aCG1CDHD2 [measure angle [list $CG1 $CD $HD2]];  
	puts   $anglemeasure "$aCG1CDHD2";
	set aCG1CDHD3 [measure angle [list $CG1 $CD $HD3]];  
	puts   $anglemeasure "$aCG1CDHD3";
	
    for {set j 0} {$j <10} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB [measure dihed [list $N $CA $CB $HB]]; 
	puts   $dihdmeasure "$dNCACBHB";
	
	set dNCACBCG2 [measure dihed [list $N $CA $CB $CG2]]; 
	puts   $dihdmeasure "$dNCACBCG2";
	
	set dCACBCG2HG21 [measure dihed [list $CA $CB $CG2 $HG21]]; 
	puts   $dihdmeasure "$dCACBCG2HG21";
	set dCACBCG2HG22 [measure dihed [list $CA $CB $CG2 $HG22]]; 
	puts   $dihdmeasure "$dCACBCG2HG22";
	set dCACBCG2HG23 [measure dihed [list $CA $CB $CG2 $HG23]]; 
	puts   $dihdmeasure "$dCACBCG2HG23";
	
	set dNCACBCG1 [measure dihed [list $N $CA $CB $CG1]]; 
	puts   $dihdmeasure "$dNCACBCG1";
	set dCACBCG1HG11 [measure dihed [list $CA $CB $CG1 $HG11]]; 
	puts   $dihdmeasure "$dCACBCG1HG11";
	set dCACBCG1HG12 [measure dihed [list $CA $CB $CG1 $HG12]]; 
	puts   $dihdmeasure "$dCACBCG1HG12";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dCACBCG1CD [measure dihed [list $CA $CB $CG1 $CD]]; 
	puts   $dihdmeasure "$dCACBCG1CD";
		
	set dCBCG1CDHD1 [measure dihed [list $CB $CG1 $CD $HD1]]; 
	puts   $dihdmeasure "$dCBCG1CDHD1";
	set dCBCG1CDHD2 [measure dihed [list $CB $CG1 $CD $HD2]]; 
	puts   $dihdmeasure "$dCBCG1CDHD2";
	set dCBCG1CDHD3 [measure dihed [list $CB $CG1 $CD $HD3]]; 
	puts   $dihdmeasure "$dCBCG1CDHD3";

	for {set j 0} {$j <10} {incr j} {
	puts   $dihdmeasure "0";
	}
    }	
	if {[string compare [lindex $res_name 0] MET]==0} {
	puts   $resi "$MET";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG2 [lindex $s_index $c_atom];
	incr c_atom;
	set SD [lindex $s_index $c_atom];
	incr c_atom;
	set CE [lindex $s_index $c_atom];
   	incr c_atom;	
	set HE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE2 [lindex $s_index $c_atom];
	incr c_atom;
	set HE3 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG1 [measure bond [list $CG $HG1]];  
	puts   $bondmeasure "$bCGHG1";
	
	set bCGHG2 [measure bond [list $CG $HG2]];  
	puts   $bondmeasure "$bCGHG2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGSD [measure bond [list $CG $SD]];  
	puts   $bondmeasure "$bCGSD";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bSDCE [measure bond [list $SD $CE]];  
	puts   $bondmeasure "$bSDCE";
	
	set bCEHE1 [measure bond [list $CE $HE1]];  
	puts   $bondmeasure "$bCEHE1";
	set bCEHE2 [measure bond [list $CE $HE2]];  
	puts   $bondmeasure "$bCEHE2";
	set bCEHE3 [measure bond [list $CE $HE3]];  
	puts   $bondmeasure "$bCEHE3";
	
	for {set j 0} {$j <7} {incr j} {
	puts   $bondmeasure "0";
    }

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG1 [measure angle [list $CB $CG $HG1]];  
	puts   $anglemeasure "$aCBCGHG1";
	
	set aCBCGHG2 [measure angle [list $CB $CG $HG2]];  
	puts   $anglemeasure "$aCBCGHG2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGSD [measure angle [list $CB $CG $SD]];  
	puts   $anglemeasure "$aCBCGSD";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
		
	set aCGSDCE [measure angle [list $CG $SD $CE]];  
	puts   $anglemeasure "$aCGSDCE";
	
	set aSDCEHE1 [measure angle [list $SD $CE $HE1]];  
	puts   $anglemeasure "$aSDCEHE1";
	
	set aSDCEHE2 [measure angle [list $SD $CE $HE2]];  
	puts   $anglemeasure "$aSDCEHE2";
	set aSDCEHE3 [measure angle [list $SD $CE $HE3]];  
	puts   $anglemeasure "$aSDCEHE3";
	
	for {set j 0} {$j <7} {incr j} {
	puts   $anglemeasure "0";
    }

	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGHG1 [measure dihed [list $CA $CB $CG $HG1]]; 
	puts   $dihdmeasure "$dCACBCGHG1";
	set dCACBCGHG2 [measure dihed [list $CA $CB $CG $HG2]]; 
	puts   $dihdmeasure "$dCACBCGHG2";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCGSD [measure dihed [list $CA $CB $CG $SD]]; 
	puts   $dihdmeasure "$dCACBCGSD";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCBCGSDCE [measure dihed [list $CB $CG $SD $CE]]; 
	puts   $dihdmeasure "$dCBCGSDCE";

	set dCGSDCEHE1 [measure dihed [list $CG $SD $CE $HE1]]; 
	puts   $dihdmeasure "$dCGSDCEHE1";
	set dCGSDCEHE2 [measure dihed [list $CG $SD $CE $HE2]]; 
	puts   $dihdmeasure "$dCGSDCEHE2";
	set dCGSDCEHE3 [measure dihed [list $CG $SD $CE $HE3]]; 
	puts   $dihdmeasure "$dCGSDCEHE3";
	
    for {set j 0} {$j <7} {incr j} {
	puts   $dihdmeasure "0";
    }

    }
	
	if {[string compare [lindex $res_name 0] GLN]==0} {
	puts   $resi "$GLN";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG2 [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];
	incr c_atom;
	set OE1 [lindex $s_index $c_atom];
   	incr c_atom;	
	set NE2 [lindex $s_index $c_atom];
	incr c_atom;
	set HE21 [lindex $s_index $c_atom];
	incr c_atom;
	set HE22 [lindex $s_index $c_atom];

	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG1 [measure bond [list $CG $HG1]];  
	puts   $bondmeasure "$bCGHG1";
	
	set bCGHG2 [measure bond [list $CG $HG2]];  
	puts   $bondmeasure "$bCGHG2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGCD [measure bond [list $CG $CD]];  
	puts   $bondmeasure "$bCGCD";
	
	set bCDOE1 [measure bond [list $CD $OE1]];  
	puts   $bondmeasure "$bCDOE1";
	puts   $bondmeasure "0";
	
	set bCDNE2 [measure bond [list $CD $NE2]];  
	puts   $bondmeasure "$bCDNE2";
	set bNE2HE21 [measure bond [list $NE2 $HE21]];  
	puts   $bondmeasure "$bNE2HE21";
	set bNE2HE22 [measure bond [list $NE2 $HE22]];  
	puts   $bondmeasure "$bNE2HE22";

	for {set j 0} {$j <8} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG1 [measure angle [list $CB $CG $HG1]];  
	puts   $anglemeasure "$aCBCGHG1";
	
	set aCBCGHG2 [measure angle [list $CB $CG $HG2]];  
	puts   $anglemeasure "$aCBCGHG2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGCD [measure angle [list $CB $CG $CD]];  
	puts   $anglemeasure "$aCBCGCD";
		
	set aCGCDOE1 [measure angle [list $CG $CD $OE1]];  
	puts   $anglemeasure "$aCGCDOE1";
	puts   $anglemeasure "0";
	
	set aCGCDNE2 [measure angle [list $CG $CD $NE2]];  
	puts   $anglemeasure "$aCGCDNE2";
	
	set aCDNE2HE21 [measure angle [list $CD $NE2 $HE21]];  
	puts   $anglemeasure "$aCDNE2HE21";
	set aCDNE2HE22 [measure angle [list $CD $NE2 $HE22]];  
	puts   $anglemeasure "$aCDNE2HE22";
	
	for {set j 0} {$j <8} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGHG1 [measure dihed [list $CA $CB $CG $HG1]]; 
	puts   $dihdmeasure "$dCACBCGHG1";
	set dCACBCGHG2 [measure dihed [list $CA $CB $CG $HG2]]; 
	puts   $dihdmeasure "$dCACBCGHG2";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCGCD [measure dihed [list $CA $CB $CG $CD]]; 
	puts   $dihdmeasure "$dCACBCGCD";
	
	set dCBCGCDOE1 [measure dihed [list $CB $CG $CD $OE1]]; 
	puts   $dihdmeasure "$dCBCGCDOE1";
	puts   $dihdmeasure "0";
	
	set dCBCGCDNE2 [measure dihed [list $CB $CG $CD $NE2]]; 
	puts   $dihdmeasure "$dCBCGCDNE2";
	
	set dCGCDNE2HE21 [measure dihed [list $CG $CD $NE2 $HE21]]; 
	puts   $dihdmeasure "$dCGCDNE2HE21";
	set dCGCDNE2HE22 [measure dihed [list $CG $CD $NE2 $HE22]]; 
	puts   $dihdmeasure "$dCGCDNE2HE22";
	
	for {set j 0} {$j <8} {incr j} {
	puts   $dihdmeasure "0";
	}	
    }
	
	if {[string compare [lindex $res_name 0] ASN]==0} {
	puts   $resi "$ASN";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set OD1 [lindex $s_index $c_atom];
	incr c_atom;
	set ND2 [lindex $s_index $c_atom];
   	incr c_atom;	
	set HD21 [lindex $s_index $c_atom];
	incr c_atom;
	set HD22 [lindex $s_index $c_atom];

	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGOD1 [measure bond [list $CG $OD1]];  
	puts   $bondmeasure "$bCGOD1";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGND2 [measure bond [list $CG $ND2]];  
	puts   $bondmeasure "$bCGND2";	
	set bND2HD21 [measure bond [list $ND2 $HD21]];  
	puts   $bondmeasure "$bND2HD21";
	set bND2HD22 [measure bond [list $ND2 $HD22]];  
	puts   $bondmeasure "$bND2HD22";

	for {set j 0} {$j <11} {incr j} {
	puts   $bondmeasure "0";
    }

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGOD1 [measure angle [list $CB $CG $OD1]];  
	puts   $anglemeasure "$aCBCGOD1";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGND2 [measure angle [list $CB $CG $ND2]];  
	puts   $anglemeasure "$aCBCGND2";
		
	set aCGND2HD21 [measure angle [list $CG $ND2 $HD21]];  
	puts   $anglemeasure "$aCGND2HD21";
	set aCGND2HD22 [measure angle [list $CG $ND2 $HD22]];  
	puts   $anglemeasure "$aCGND2HD22";
	
	for {set j 0} {$j <11} {incr j} {
	puts   $anglemeasure "0";
    }
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGOD1 [measure dihed [list $CA $CB $CG $OD1]]; 
	puts   $dihdmeasure "$dCACBCGOD1";
	
	puts   $dihdmeasure "0";	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	
	set dCACBCGND2 [measure dihed [list $CA $CB $CG $ND2]]; 
	puts   $dihdmeasure "$dCACBCGND2";
	set dCBCGND2HG21 [measure dihed [list $CB $CG $ND2 $HG21]]; 
	puts   $dihdmeasure "$dCBCGND2HG21";
	set dCBCGND2HG22 [measure dihed [list $CB $CG $ND2 $HG22]]; 
	puts   $dihdmeasure "$dCBCGND2HG22";

	for {set j 0} {$j <11} {incr j} {
	puts   $dihdmeasure "0";
    }

    }
	
	if {[string compare [lindex $res_name 0] VAL]==0} {
	puts   $resi "$VAL";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB [lindex $s_index $c_atom];
	incr c_atom;
	set CG1 [lindex $s_index $c_atom];
	incr c_atom;
	set CG2 [lindex $s_index $c_atom];
	incr c_atom;
	incr c_atom;
	incr c_atom;
	set HG11 [lindex $s_index $c_atom];
	incr c_atom;
	set HG12 [lindex $s_index $c_atom];
	incr c_atom;
	set HG13 [lindex $s_index $c_atom];
	incr c_atom;
	set HG21 [lindex $s_index $c_atom];
	incr c_atom;
	set HG22 [lindex $s_index $c_atom];
	incr c_atom;
	set HG23 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB [measure bond [list $CB $HB]];   
	puts   $bondmeasure "$bCBHB";
	set bCBCG1 [measure bond [list $CB $CG1]];   
	puts   $bondmeasure "$bCBCG1";
	
	set bCG1HG11 [measure bond [list $CG1 $HG11]];   
	puts   $bondmeasure "$bCG1HG11";
	
	set bCG1HG12 [measure bond [list $CG1 $HG12]];   
	puts   $bondmeasure "$bCG1HG12";
	
	set bCG1HG13 [measure bond [list $CG1 $HG13]];   
	puts   $bondmeasure "$bCG1HG13";
	
	set bCBCG2 [measure bond [list $CB $CG2]];   
	puts   $bondmeasure "$bCBCG2";
	
	set bCG2HG21 [measure bond [list $CG2 $HG21]];   
	puts   $bondmeasure "$bCG2HG21";
	
	set bCG2HG22 [measure bond [list $CG2 $HG22]];   
	puts   $bondmeasure "$bCG2HG22";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	set bCG2CG23 [measure bond [list $CG2 $HG23]];   
	puts   $bondmeasure "$bCG2HG23";

	for {set j 0} {$j <13} {incr j} {
	puts   $bondmeasure "0";	
	}
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB [measure angle [list $CA $CB $HB]];  
	puts   $anglemeasure "$aCACBHB";	
	
	set aCACBCG1 [measure angle [list $CA $CB $CG1]];  
	puts   $anglemeasure "$aCACBCG1";
	set aCBCG1HG11 [measure angle [list $CB $CG1 $HG11]];  
	puts   $anglemeasure "$aCBCG1HG11";
	set aCBCG1HG12 [measure angle [list $CB $CG1 $HG12]];  
	puts   $anglemeasure "$aCBCG1HG12";
	set aCBCG1HG13 [measure angle [list $CB $CG1 $HG13]];  
	puts   $anglemeasure "$aCBCG1HG13";
	
	set aCACBCG2 [measure angle [list $CA $CB $CG2]];  
	puts   $anglemeasure "$aCACBCG2";
	set aCBCG2HG21 [measure angle [list $CB $CG2 $HG21]];  
	puts   $anglemeasure "$aCBCG2HG21";
	set aCBCG2HG22 [measure angle [list $CB $CG2 $HG22]];  
	puts   $anglemeasure "$aCBCG2HG22";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCG2HG23 [measure angle [list $CB $CG2 $HG23]];  
	puts   $anglemeasure "$aCBCG2HG23";
	
    for {set j 0} {$j <13} {incr j} {
	puts   $anglemeasure "0";
    }
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
	set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB [measure dihed [list $N $CA $CB $HB]]; 
	puts   $dihdmeasure "$dNCACBHB";

	set dNCACBCG1 [measure dihed [list $N $CA $CB $CG1]]; 
	puts   $dihdmeasure "$dNCACBCG1";
	set dCACBCG1HG11 [measure dihed [list $CA $CB $CG1 $HG11]]; 
	puts   $dihdmeasure "$dCACBCG1HG11";
	set dCACBCG1HG12 [measure dihed [list $CA $CB $CG1 $HG12]]; 
	puts   $dihdmeasure "$dCACBCG1HG12";
	set dCACBCG1HG13 [measure dihed [list $CA $CB $CG1 $HG13]]; 
	puts   $dihdmeasure "$dCACBCG1HG13";
	
	set dNCACBCG2 [measure dihed [list $N $CA $CB $CG2]]; 
	puts   $dihdmeasure "$dNCACBCG2";
	
	set dCACBCG2HG21 [measure dihed [list $CA $CB $CG2 $HG21]]; 
	puts   $dihdmeasure "$dCACBCG2HG21";
	set dCACBCG2HG22 [measure dihed [list $CA $CB $CG2 $HG22]]; 
	puts   $dihdmeasure "$dCACBCG2HG22";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	set dCACBCG2HG23 [measure dihed [list $CA $CB $CG2 $HG23]]; 
	puts   $dihdmeasure "$dCACBCG2HG23";
	
	for {set j 0} {$j <13} {incr j} {
	puts   $dihdmeasure "0";
	}
    }	
	if {[string compare [lindex $res_name 0] ARG]==0} {
	puts   $resi "$ARG";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set HG1 [lindex $s_index $c_atom];
	incr c_atom;
	set HG2 [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
   	incr c_atom;
	set NE [lindex $s_index $c_atom];
   	incr c_atom;	
	set HE [lindex $s_index $c_atom];
	incr c_atom;
	set CZ [lindex $s_index $c_atom];
	incr c_atom;
	set NH1 [lindex $s_index $c_atom];
	incr c_atom;
	set HH11 [lindex $s_index $c_atom];
	incr c_atom;
	set HH12 [lindex $s_index $c_atom];
	incr c_atom;
	set NH2 [lindex $s_index $c_atom];
	incr c_atom;
	set HH21 [lindex $s_index $c_atom];
	incr c_atom;
	set HH22 [lindex $s_index $c_atom];
	
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	
	set bCGHG1 [measure bond [list $CG $HG1]];  
	puts   $bondmeasure "$bCGHG1";
	
	set bCGHG2 [measure bond [list $CG $HG2]];  
	puts   $bondmeasure "$bCGHG2";
	
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	puts   $bondmeasure "0";
	
	set bCGCD [measure bond [list $CG $CD]];  
	puts   $bondmeasure "$bCGCD";
	
	set bCDHD1 [measure bond [list $CD $HD1]];  
	puts   $bondmeasure "$bCDHD1";
	set bCDHD2 [measure bond [list $CD $HD2]];  
	puts   $bondmeasure "$bCDHD2";
	
	set bCDNE [measure bond [list $CD $NE]];  
	puts   $bondmeasure "$bCDNE";
	
	set bNEHE [measure bond [list $NE $HE]];  
	puts   $bondmeasure "$bNEHE";
	puts   $bondmeasure "0";
	
	set bNECZ [measure bond [list $NE $CZ]];  
	puts   $bondmeasure "$bNECZ";
	set bCZNH1 [measure bond [list $CZ $NH1]];  
	puts   $bondmeasure "$bCZNH1";
	puts   $bondmeasure "0";
	
	set bNH1HH11 [measure bond [list $NH1 $HH11]];  
	puts   $bondmeasure "$bNH1HH11";
	set bNH1HH12 [measure bond [list $NH1 $HH12]];  
	puts   $bondmeasure "$bNH1HH12";

	set bCZNH2 [measure bond [list $CZ $NH2]];  
	puts   $bondmeasure "$bCZNH2";	
	set bNH2HH21 [measure bond [list $NH2 $HH21]];  
	puts   $bondmeasure "$bNH2HH21";
	set bNH2HH22 [measure bond [list $NH2 $HH22]];  
	puts   $bondmeasure "$bNH2HH22";
	
	
	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	
    puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	
	set aCBCGHG1 [measure angle [list $CB $CG $HG1]];  
	puts   $anglemeasure "$aCBCGHG1";
	
	set aCBCGHG2 [measure angle [list $CB $CG $HG2]];  
	puts   $anglemeasure "$aCBCGHG2";
	
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	puts   $anglemeasure "0";
	
	set aCBCGCD [measure angle [list $CB $CG $CD]];  
	puts   $anglemeasure "$aCBCGCD";
		
	set aCGCDHD1 [measure angle [list $CG $CD $HD1]];  
	puts   $anglemeasure "$aCGCDHD1";
	
	set aCGCDHD2 [measure angle [list $CG $CD $HD2]];  
	puts   $anglemeasure "$aCGCDHD2";
	
	set aCGCDNE [measure angle [list $CG $CD $NE]];  
	puts   $anglemeasure "$aCGCDNE";
	
	set aCDNEHE [measure angle [list $CD $NE $HE]];  
	puts   $anglemeasure "$aCDNEHE";
	puts   $anglemeasure "0";
	
	set aCDNECZ [measure angle [list $CD $NE $CZ]];  
	puts   $anglemeasure "$aCDNECZ";
	
	set aNECZNH1 [measure angle [list $NE $CZ $NH1]];  
	puts   $anglemeasure "$aNECZNH1";
	puts   $anglemeasure "0";
	set aCZNH1HH11 [measure angle [list $CZ $NH1 $HH11]];  
	puts   $anglemeasure "$aCZNH1HH11";
	set aCZNH1HH12 [measure angle [list $CZ $NH1 $HH12]];  
	puts   $anglemeasure "$aCZNH1HH12";
	set aNECZNH2 [measure angle [list $NE $CZ $NH2]];  
	puts   $anglemeasure "$aNECZNH2";
	set aCZNH2HH21 [measure angle [list $CZ $NH2 $HH21]];  
	puts   $anglemeasure "$aCZNH2HH21";
	set aCZNH2HH22 [measure angle [list $CZ $NH2 $HH22]];  
	puts   $anglemeasure "$aCZNH2HH22";

	

	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";

	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
		
	set dCACBCGHG1 [measure dihed [list $CA $CB $CG $HG1]]; 
	puts   $dihdmeasure "$dCACBCGHG1";
	
	set dCACBCGHG2 [measure dihed [list $CA $CB $CG $HG2]]; 
	puts   $dihdmeasure "$dCACBCGHG2";
	
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	puts   $dihdmeasure "0";
	set dCACBCGCD [measure dihed [list $CA $CB $CG $CD]]; 
	puts   $dihdmeasure "$dCACBCGCD";
	
	set dCBCGCDHD1 [measure dihed [list $CB $CG $CD $HD1]]; 
	puts   $dihdmeasure "$dCBCGCDHD1";
	
	set dCBCGCDHD2 [measure dihed [list $CB $CG $CD $HD2]]; 
	puts   $dihdmeasure "$dCBCGCDHD2";
	
	set dCBCGCDNE [measure dihed [list $CB $CG $CD $NE]]; 
	puts   $dihdmeasure "$dCBCGCDNE";
	
	set dCGCDNEHE [measure dihed [list $CG $CD $NE $HE]]; 
	puts   $dihdmeasure "$dCGCDNEHE";
	puts   $dihdmeasure "0";
	
	set dCGCDNECZ [measure dihed [list $CG $CD $NE $CZ]]; 
	puts   $dihdmeasure "$dCGCDNECZ";
	
	set dCDNECZNH1 [measure dihed [list $CD $NE $CZ $NH1]]; 
	puts   $dihdmeasure "$dCDNECZNH1";
	puts   $dihdmeasure "0";
	set dNECZNH1HH11 [measure dihed [list $NE $CZ $NH1 $HH11]]; 
	puts   $dihdmeasure "$dNECZNH1HH11";
	set dNECZNH1HH12 [measure dihed [list $NE $CZ $NH1 $HH12]]; 
	puts   $dihdmeasure "$dNECZNH1HH12";
	
	set dCDNECZNH2 [measure dihed [list $CD $NE $CZ $NH2]]; 
	puts   $dihdmeasure "$dCDNECZNH2";
	set dNECZNH2HH21 [measure dihed [list $NE $CZ $NH2 $HH21]]; 
	puts   $dihdmeasure "$dNECZNH2HH21";
	set dNECZNH2HH22 [measure dihed [list $NE $CZ $NH2 $HH22]]; 
	puts   $dihdmeasure "$dNECZNH2HH22";

}
	if {[string compare [lindex $res_name 0] PRO]==0} {
	puts   $resi "$PRO";
	
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD [lindex $s_index $c_atom];


	
	puts [$side_i get {name index}]
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

   set bCACB [measure bond [list $CA $CB]];
	puts  $bondmeasure "$bCACB"
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGCD [measure bond [list $CG $CD]];  
	puts   $bondmeasure "$bCGCD";
	
	for {set j 0} {$j <23} {incr j} {
	puts   $bondmeasure "0";		
	}
	

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGCD [measure angle [list $CB $CG $CD]];  
	puts   $anglemeasure "$aCBCGCD";
	
	for {set j 0} {$j <23} {incr j} {
	puts   $anglemeasure "0";		
	}		
    set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";	
	set dCACBCGCD [measure dihed [list $CA $CB $CG $CD]]; 
	puts   $dihdmeasure "$dCACBCGCD";
	
	for {set j 0} {$j <23} {incr j} {
	puts   $dihdmeasure "0";	
	}
}
	if {[string compare [lindex $res_name 0] TYR]==0} {
	puts   $resi "$TYR";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set CE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set CZ [lindex $s_index $c_atom];
   	incr c_atom;
	set OH [lindex $s_index $c_atom];
   	incr c_atom;	
	set HH [lindex $s_index $c_atom];
	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
	incr c_atom;
	set CE2 [lindex $s_index $c_atom];
	incr c_atom;
	set HE2 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB]];
	puts  $bondmeasure "$bCACB"
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGCD1 [measure bond [list $CG $CD1]];  
	puts   $bondmeasure "$bCGCD1";
	set bCD1HD1 [measure bond [list $CD1 $HD1]];  
	puts   $bondmeasure "$bCD1HD1";
	set bCD1CE1 [measure bond [list $CD1 $CE1]];  
	puts   $bondmeasure "$bCD1CE1";
	set bCE1HE1 [measure bond [list $CE1 $HE1]];  
	puts   $bondmeasure "$bCE1HE1";
	set bCE1CZ [measure bond [list $CE1 $CZ]];  
	puts   $bondmeasure "$bCE1CZ";
	set bCZOH [measure bond [list $CZ $OH]];  
	puts   $bondmeasure "$bCZOH";
	set bOHHH [measure bond [list $OH $HH]];  
	puts   $bondmeasure "$bOHHH";
	set bCZCE2 [measure bond [list $CZ $CE2]];  
	puts   $bondmeasure "$bCZCE2";
	set bCE2HE2 [measure bond [list $CE2 $HE2]];  
	puts   $bondmeasure "$bCE2HE2";
	set bCE2CD2 [measure bond [list $CE2 $CD2]];  
	puts   $bondmeasure "$bCE2CD2";
	set bCD2HD2 [measure bond [list $CD2 $HD2]];  
	puts   $bondmeasure "$bCD2HD2";

	
	for {set j 0} {$j <11} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGCD1 [measure angle [list $CB $CG $CD1]];  
	puts   $anglemeasure "$aCBCGCD1";
	set aCGCD1HD1 [measure angle [list $CG $CD1 $HD1]];  
	puts   $anglemeasure "$aCGCD1HD1";
	set aCGCD1CE1 [measure angle [list $CG $CD1 $CE1]];  
	puts   $anglemeasure "$aCGCD1CE1";
	set aCD1CE1HE1 [measure angle [list $CD1 $CE1 $HE1]];  
	puts   $anglemeasure "$aCD1CE1HE1";
	set aCD1CE1CZ [measure angle [list $CD1 $CE1 $CZ]];  
	puts   $anglemeasure "$aCD1CE1CZ";
	set aCE1CZOH [measure angle [list $CE1 $CZ $OH]];  
	puts   $anglemeasure "$aCE1CZOH";
	set aCZOHHH [measure angle [list $CZ $OH $HH]];  
	puts   $anglemeasure "$aCZOHHH";
	set aCE1CZCE2 [measure angle [list $CE1 $CZ $CE2]];  
	puts   $anglemeasure "$aCE1CZCE2";
	set aCZCE2HE2 [measure angle [list $CZ $CE2 $HE2]];  
	puts   $anglemeasure "$aCZCE2HE2";
	set aCZCE2CD2 [measure angle [list $CZ $CE2 $CD2]];  
	puts   $anglemeasure "$aCZCE2CD2";
	set aCE2CD2HD2 [measure angle [list $CE2 $CD2 $HD2]];  
	puts   $anglemeasure "$aCE2CD2HD2";
		
	for {set j 0} {$j <11} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
	
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	set dCACBCGCD1 [measure dihed [list $CA $CB $CG $CD1]]; 
	puts   $dihdmeasure "$dCACBCGCD1";
	set dCBCGCD1HD1 [measure dihed [list $CB $CG $CD1 $HD1]]; 
	puts   $dihdmeasure "$dCBCGCD1HD1";
	set dCBCGCD1CE1 [measure dihed [list $CB $CG $CD1 $CE1]]; 
	puts   $dihdmeasure "$dCBCGCD1CE1";
	set dCGCD1CE1HE1 [measure dihed [list $CG $CD1 $CE1 $HE1]]; 
	puts   $dihdmeasure "$dCGCD1CE1HE1";
	set dCGCD1CE1CZ [measure dihed [list $CG $CD1 $CE1 $CZ]]; 
	puts   $dihdmeasure "$dCGCD1CE1CZ";
	set dCD1CE1CZOH [measure dihed [list $CD1 $CE1 $CZ $OH]]; 
	puts   $dihdmeasure "$dCD1CE1CZOH";
	set dCE1CZOHHH [measure dihed [list $CE1 $CZ $OH $HH]]; 
	puts   $dihdmeasure "$dCE1CZOHHH";
	set dCD1CE1CZCE2 [measure dihed [list $CD1 $CE1 $CZ $CE2]]; 
	puts   $dihdmeasure "$dCD1CE1CZCE2";
	set dCE1CZCE2HE2 [measure dihed [list $CE1 $CZ $CE2 $HE2]]; 
	puts   $dihdmeasure "$dCE1CZCE2HE2";
	set dCE1CZCE2CD2 [measure dihed [list $CE1 $CZ $CE2 $CD2]]; 
	puts   $dihdmeasure "$dCE1CZCE2CD2";
	set dCZCE2CD2HD2 [measure dihed [list $CZ $CE2 $CD2 $HD2]]; 
	puts   $dihdmeasure "$dCZCE2CD2HD2";
	
	for {set j 0} {$j <11} {incr j} {	
	puts   $dihdmeasure "0";
    }
}
	if {[string compare [lindex $res_name 0] PHE]==0} {
	puts   $resi "$PHE";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set CE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set CZ [lindex $s_index $c_atom];
   	incr c_atom;
	set HZ [lindex $s_index $c_atom];
   	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
	incr c_atom;
	set CE2 [lindex $s_index $c_atom];
	incr c_atom;
	set HE2 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGCD1 [measure bond [list $CG $CD1]];  
	puts   $bondmeasure "$bCGCD1";
	set bCD1HD1 [measure bond [list $CD1 $HD1]];  
	puts   $bondmeasure "$bCD1HD1";
	set bCD1CE1 [measure bond [list $CD1 $CE1]];  
	puts   $bondmeasure "$bCD1CE1";
	set bCE1HE1 [measure bond [list $CE1 $HE1]];  
	puts   $bondmeasure "$bCE1HE1";
	set bCE1CZ [measure bond [list $CE1 $CZ]];  
	puts   $bondmeasure "$bCE1CZ";
	set bCZHZ [measure bond [list $CZ $HZ]];  
	puts   $bondmeasure "$bCZHZ";
	set bCZCE2 [measure bond [list $CZ $CE2]];  
	puts   $bondmeasure "$bCZCE2";
	set bCE2HE2 [measure bond [list $CE2 $HE2]];  
	puts   $bondmeasure "$bCE2HE2";
	set bCE2CD2 [measure bond [list $CE2 $CD2]];  
	puts   $bondmeasure "$bCE2CD2";
	set bCD2HD2 [measure bond [list $CD2 $HD2]];  
	puts   $bondmeasure "$bCD2HD2";

	
	for {set j 0} {$j <12} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGCD1 [measure angle [list $CB $CG $CD1]];  
	puts   $anglemeasure "$aCBCGCD1";
	set aCGCD1HD1 [measure angle [list $CG $CD1 $HD1]];  
	puts   $anglemeasure "$aCGCD1HD1";
	set aCGCD1CE1 [measure angle [list $CG $CD1 $CE1]];  
	puts   $anglemeasure "$aCGCD1CE1";
	set aCD1CE1HE1 [measure angle [list $CD1 $CE1 $HE1]];  
	puts   $anglemeasure "$aCD1CE1HE1";
	set aCD1CE1CZ [measure angle [list $CD1 $CE1 $CZ]];  
	puts   $anglemeasure "$aCD1CE1CZ";
	set aCE1CZHZ [measure angle [list $CE1 $CZ $HZ]];  
	puts   $anglemeasure "$aCE1CZHZ";
	set aCE1CZCE2 [measure angle [list $CE1 $CZ $CE2]];  
	puts   $anglemeasure "$aCE1CZCE2";
	set aCZCE2HE2 [measure angle [list $CZ $CE2 $HE2]];  
	puts   $anglemeasure "$aCZCE2HE2";
	set aCZCE2CD2 [measure angle [list $CZ $CE2 $CD2]];  
	puts   $anglemeasure "$aCZCE2CD2";
	set aCE2CD2HD2 [measure angle [list $CE2 $CD2 $HD2]];  
	puts   $anglemeasure "$aCE2CD2HD2";
		
	for {set j 0} {$j <12} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	set dCACBCGCD1 [measure dihed [list $CA $CB $CG $CD1]]; 
	puts   $dihdmeasure "$dCACBCGCD1";
	set dCBCGCD1HD1 [measure dihed [list $CB $CG $CD1 $HD1]]; 
	puts   $dihdmeasure "$dCBCGCD1HD1";
	set dCBCGCD1CE1 [measure dihed [list $CB $CG $CD1 $CE1]]; 
	puts   $dihdmeasure "$dCBCGCD1CE1";
	set dCGCD1CE1HE1 [measure dihed [list $CG $CD1 $CE1 $HE1]]; 
	puts   $dihdmeasure "$dCGCD1CE1HE1";
	set dCGCD1CE1CZ [measure dihed [list $CG $CD1 $CE1 $CZ]]; 
	puts   $dihdmeasure "$dCGCD1CE1CZ";
	set dCD1CE1CZHZ [measure dihed [list $CD1 $CE1 $CZ $HZ]]; 
	puts   $dihdmeasure "$dCD1CE1CZHZ";
	set dCD1CE1CZCE2 [measure dihed [list $CD1 $CE1 $CZ $CE2]]; 
	puts   $dihdmeasure "$dCD1CE1CZCE2";
	set dCE1CZCE2HE2 [measure dihed [list $CE1 $CZ $CE2 $HE2]]; 
	puts   $dihdmeasure "$dCE1CZCE2HE2";
	set dCE1CZCE2CD2 [measure dihed [list $CE1 $CZ $CE2 $CD2]]; 
	puts   $dihdmeasure "$dCE1CZCE2CD2";
	set dCZCE2CD2HD2 [measure dihed [list $CZ $CE2 $CD2 $HD2]]; 
	puts   $dihdmeasure "$dCZCE2CD2HD2";
	
	for {set j 0} {$j <12} {incr j} {	
	puts   $dihdmeasure "0";
    }
}

	if {[string compare [lindex $res_name 0] TRP]==0} {
	puts   $resi "$TRP";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CD1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set NE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set CE2 [lindex $s_index $c_atom];
   	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
   	incr c_atom;
	set CE3 [lindex $s_index $c_atom];
   	incr c_atom;
	set HE3 [lindex $s_index $c_atom];
   	incr c_atom;
	set CZ3 [lindex $s_index $c_atom];
   	incr c_atom;
	set HZ3 [lindex $s_index $c_atom];
   	incr c_atom;
	set CZ2 [lindex $s_index $c_atom];
	incr c_atom;
	set HZ2 [lindex $s_index $c_atom];
	incr c_atom;
	set CH2 [lindex $s_index $c_atom];
	incr c_atom;
	set HH2 [lindex $s_index $c_atom];
	
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGCD1 [measure bond [list $CG $CD1]];  
	puts   $bondmeasure "$bCGCD1";
	set bCD1HD1 [measure bond [list $CD1 $HD1]];  
	puts   $bondmeasure "$bCD1HD1";
	set bCD1NE1 [measure bond [list $CD1 $NE1]];  
	puts   $bondmeasure "$bCD1NE1";
	set bNE1HE1 [measure bond [list $NE1 $HE1]];  
	puts   $bondmeasure "$bNE1HE1";
	set bNE1CE2 [measure bond [list $NE1 $CE2]];  
	puts   $bondmeasure "$bNE1CE2";
	set bCE2CD2 [measure bond [list $CE2 $CD2]];  
	puts   $bondmeasure "$bCE2CD2";
	set bCD2CE3 [measure bond [list $CD2 $CE3]];  
	puts   $bondmeasure "$bCD2CE3";
	set bCE3HE3 [measure bond [list $CE3 $HE3]];  
	puts   $bondmeasure "$bCE3HE3";
	set bCE3CZ3 [measure bond [list $CE3 $CZ3]];  
	puts   $bondmeasure "$bCE3CZ3";
	set bCZ3HZ3 [measure bond [list $CZ3 $HZ3]];  
	puts   $bondmeasure "$bCZ3HZ3";
	set bCZ3CH2 [measure bond [list $CZ3 $CH2]];  
	puts   $bondmeasure "$bCZ3CH2";
	set bCH2HH2 [measure bond [list $CH2 $HH2]];  
	puts   $bondmeasure "$bCH2HH2";
	set bCH2CZ2 [measure bond [list $CH2 $CZ2]];  
	puts   $bondmeasure "$bCH2CZ2";
	set bCZ2HZ2 [measure bond [list $CZ2 $HZ2]];  
	puts   $bondmeasure "$bCZ2HZ2";

	
	for {set j 0} {$j <8} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGCD1 [measure angle [list $CB $CG $CD1]];  
	puts   $anglemeasure "$aCBCGCD1";
	set aCGCD1HD1 [measure angle [list $CG $CD1 $HD1]];  
	puts   $anglemeasure "$aCGCD1HD1";
	set aCGCD1NE1 [measure angle [list $CG $CD1 $NE1]];  
	puts   $anglemeasure "$aCGCD1NE1";
	set aCD1NE1HE1 [measure angle [list $CD1 $NE1 $HE1]];  
	puts   $anglemeasure "$aCD1NE1HE1";
	set aCD1NE1CE2 [measure angle [list $CD1 $NE1 $CE2]];  
	puts   $anglemeasure "$aCD1NE1CE2";
	set aNE1CE2CD2 [measure angle [list $NE1 $CE2 $CD2]];  
	puts   $anglemeasure "$aNE1CE2CD2";
	set aCE2CD2CE3 [measure angle [list $CE2 $CD2 $CE3]];  
	puts   $anglemeasure "$aCE2CD2CE3";
	set aCD2CE3HE3 [measure angle [list $CD2 $CE3 $HE3]];  
	puts   $anglemeasure "$aCD2CE3HE3";
	set aCD2CE3CZ3 [measure angle [list $CD2 $CE3 $CZ3]];  
	puts   $anglemeasure "$aCD2CE3CZ3";
	set aCE3CZ3HZ3 [measure angle [list $CE3 $CZ3 $HZ3]];  
	puts   $anglemeasure "$aCE3CZ3HZ3";
	set aCE3CZ3CH2 [measure angle [list $CE3 $CZ3 $CH2]];  
	puts   $anglemeasure "$aCE3CZ3CH2";
	set aCZ3CH2HH2 [measure angle [list $CZ3 $CH2 $HH2]];  
	puts   $anglemeasure "$aCZ3CH2HH2";
	set aCZ3CH2CZ2 [measure angle [list $CZ3 $CH2 $CZ2]];  
	puts   $anglemeasure "$aCZ3CH2CZ2";
	set aCH2CZ2HZ2 [measure angle [list $CH2 $CZ2 $HZ2]];  
	puts   $anglemeasure "$aCH2CZ2HZ2";	
		
	for {set j 0} {$j <8} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	set dCACBCGCD1 [measure dihed [list $CA $CB $CG $CD1]]; 
	puts   $dihdmeasure "$dCACBCGCD1";
	set dCBCGCD1HD1 [measure dihed [list $CB $CG $CD1 $HD1]]; 
	puts   $dihdmeasure "$dCBCGCD1HD1";
	set dCBCGCD1NE1 [measure dihed [list $CB $CG $CD1 $NE1]]; 
	puts   $dihdmeasure "$dCBCGCD1NE1";
	set dCGCD1NE1HE1 [measure dihed [list $CG $CD1 $NE1 $HE1]]; 
	puts   $dihdmeasure "$dCGCD1NE1HE1";
	set dCGCD1NE1CE2 [measure dihed [list $CG $CD1 $NE1 $CE2]]; 
	puts   $dihdmeasure "$dCGCD1NE1CE2";
	set dCD1NE1CE2CD2 [measure dihed [list $CD1 $NE1 $CE2 $CD2]]; 
	puts   $dihdmeasure "$dCD1NE1CE2CD2";
	set dNE1CE2CD2CE3 [measure dihed [list $NE1 $CE2 $CD2 $CE3]]; 
	puts   $dihdmeasure "$dNE1CE2CD2CE3";
	set dCE2CD2CE3HE3 [measure dihed [list $CE2 $CD2 $CE3 $HE3]]; 
	puts   $dihdmeasure "$dCE2CD2CE3HE3";
	set dCE2CD2CE3CZ3 [measure dihed [list $CE2 $CD2 $CE3 $CZ3]]; 
	puts   $dihdmeasure "$dCE2CD2CE3CZ3";
	set dCD2CE3CZ3HZ3 [measure dihed [list $CD2 $CE3 $CZ3 $HZ3]]; 
	puts   $dihdmeasure "$dCD2CE3CZ3HZ3";
	set dCD2CE3CZ3CH2 [measure dihed [list $CD2 $CE3 $CZ3 $CH2]]; 
	puts   $dihdmeasure "$dCD2CE3CZ3CH2";
	set dCE3CZ3CH2HH2 [measure dihed [list $CE3 $CZ3 $CH2 $HH2]]; 
	puts   $dihdmeasure "$dCE3CZ3CH2HH2";
	set dCE3CZ3CH2CZ2 [measure dihed [list $CE3 $CZ3 $CH2 $CZ2]]; 
	puts   $dihdmeasure "$dCE3CZ3CH2CZ2";
	set dCZ3CH2CZ2HZ2 [measure dihed [list $CZ3 $CH2 $CZ2 $HZ2]]; 
	puts   $dihdmeasure "$dCZ3CH2CZ2HZ2";
	
	for {set j 0} {$j <8} {incr j} {	
	puts   $dihdmeasure "0";
    }
}

	if {[string compare [lindex $res_name 0] HSD]==0} {
	puts   $resi "$HSD";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set ND1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set NE2 [lindex $s_index $c_atom];
   	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGND1 [measure bond [list $CG $ND1]];  
	puts   $bondmeasure "$bCGND1";
	set bND1HD1 [measure bond [list $ND1 $HD1]];  
	puts   $bondmeasure "$bND1HD1";
	set bND1CE1 [measure bond [list $ND1 $CE1]];  
	puts   $bondmeasure "$bND1CE1";
	set bCE1HE1 [measure bond [list $NE1 $HE1]];  
	puts   $bondmeasure "$bNE1HE1";
	set bCE1NE2 [measure bond [list $CE1 $NE2]];  
	puts   $bondmeasure "$bCE1NE2";
	set bNE2CD2 [measure bond [list $NE2 $CD2]];  
	puts   $bondmeasure "$bNE2CD2";
	set bCD2HD2 [measure bond [list $CD2 $HD2]];  
	puts   $bondmeasure "$bCD2HD2";

	for {set j 0} {$j <15} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGND1 [measure angle [list $CB $CG $ND1]];  
	puts   $anglemeasure "$aCBCGND1";
	set aCGND1HD1 [measure angle [list $CG $ND1 $HD1]];  
	puts   $anglemeasure "$aCGND1HD1";
	set aCGND1CE1 [measure angle [list $CG $ND1 $CE1]];  
	puts   $anglemeasure "$aCGND1CE1";
	set aND1CE1HE1 [measure angle [list $ND1 $CE1 $HE1]];  
	puts   $anglemeasure "$aND1CE1HE1";
	set aND1CE1NE2 [measure angle [list $ND1 $CE1 $NE2]];  
	puts   $anglemeasure "$aND1CE1NE2";
	set aCE1NE2CD2 [measure angle [list $CE1 $NE2 $CD2]];  
	puts   $anglemeasure "$aCE1NE2CD2";
	set aNE2CD2HD2 [measure angle [list $NE2 $CD2 $HD2]];  
	puts   $anglemeasure "$aNE2CD2HD2";

	for {set j 0} {$j <15} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	set dCACBCGND1 [measure dihed [list $CA $CB $CG $ND1]]; 
	puts   $dihdmeasure "$dCACBCGND1";
	set dCBCGND1HD1 [measure dihed [list $CB $CG $ND1 $HD1]]; 
	puts   $dihdmeasure "$dCBCGND1HD1";
	set dCBCGND1CE1 [measure dihed [list $CB $CG $ND1 $CE1]]; 
	puts   $dihdmeasure "$dCBCGND1CE1";
	set dCGND1CE1HE1 [measure dihed [list $CG $ND1 $CE1 $HE1]]; 
	puts   $dihdmeasure "$dCGND1CE1HE1";
	set dCGND1CE1NE2 [measure dihed [list $CG $ND1 $CE1 $NE2]]; 
	puts   $dihdmeasure "$dCGND1CE1NE2";
	set dND1CE1NE2CD2 [measure dihed [list $ND1 $CE1 $NE2 $CD2]]; 
	puts   $dihdmeasure "$dND1CE1NE2CD2";
	set dCE1NE2CD2HD2 [measure dihed [list $CE1 $NE2 $CD2 $HD2]]; 
	puts   $dihdmeasure "$dCE1NE2CD2HD2";
	
	for {set j 0} {$j <15} {incr j} {	
	puts   $dihdmeasure "0";
    }
}
	if {[string compare [lindex $res_name 0] HIS]==0} {
	puts   $resi "$HIS";
	set c_atom 4
	set CB [lindex $s_index $c_atom];
	incr c_atom;
	set HB1 [lindex $s_index $c_atom];
	incr c_atom;
	set HB2 [lindex $s_index $c_atom];
	incr c_atom;
	set ND1 [lindex $s_index $c_atom];
	incr c_atom;
	set HD1 [lindex $s_index $c_atom];
	incr c_atom;
	set CG [lindex $s_index $c_atom];
	incr c_atom;
	set CE1 [lindex $s_index $c_atom];
	incr c_atom;
	set HE1 [lindex $s_index $c_atom];
   	incr c_atom;
	set NE2 [lindex $s_index $c_atom];
   	incr c_atom;
	set CD2 [lindex $s_index $c_atom];
   	incr c_atom;
	set HD2 [lindex $s_index $c_atom];
	
 	puts [$side_i get {name index}]
	
	set CA  [lindex $selatomCA [expr {$i-1}]]
	set N   [lindex $selatomN [expr {$i-1}]]
	set C   [lindex $selatomC [expr {$i-1}]]

	
    set bCACB [measure bond [ list $CA $CB ]];
	puts  $bondmeasure "$bCACB"
	set bCBHB1 [measure bond [list $CB $HB1]];   
	puts   $bondmeasure "$bCBHB1";
	set bCBHB2 [measure bond [list $CB $HB2]];   
	puts   $bondmeasure "$bCBHB2";
    set bCBCG [measure bond [list $CB $CG]];  
	puts   $bondmeasure "$bCBCG";
	set bCGND1 [measure bond [list $CG $ND1]];  
	puts   $bondmeasure "$bCGND1";
	set bND1HD1 [measure bond [list $ND1 $HD1]];  
	puts   $bondmeasure "$bND1HD1";
	set bND1CE1 [measure bond [list $ND1 $CE1]];  
	puts   $bondmeasure "$bND1CE1";
	set bCE1HE1 [measure bond [list $CE1 $HE1]];  
	puts   $bondmeasure "$bCE1HE1";
	set bCE1NE2 [measure bond [list $CE1 $NE2]];  
	puts   $bondmeasure "$bCE1NE2";
	set bNE2CD2 [measure bond [list $NE2 $CD2]];  
	puts   $bondmeasure "$bNE2CD2";
	set bCD2HD2 [measure bond [list $CD2 $HD2]];  
	puts   $bondmeasure "$bCD2HD2";

	for {set j 0} {$j <15} {incr j} {
	puts   $bondmeasure "0";
	}

	set aNCACB [measure angle [list $N $CA $CB]];  
	puts   $anglemeasure "$aNCACB";
	set aCACBHB1 [measure angle [list $CA $CB $HB1]];  
	puts   $anglemeasure "$aCACBHB1";	
	set aCACBHB2 [measure angle [list $CA $CB $HB2]];  
	puts   $anglemeasure "$aCACBHB2";
	set aCACBCG [measure angle [list $CA $CB $CG]];   
	puts   $anglemeasure "$aCACBCG";
	set aCBCGND1 [measure angle [list $CB $CG $ND1]];  
	puts   $anglemeasure "$aCBCGND1";
	set aCGND1HD1 [measure angle [list $CG $ND1 $HD1]];  
	puts   $anglemeasure "$aCGND1HD1";
	set aCGND1CE1 [measure angle [list $CG $ND1 $CE1]];  
	puts   $anglemeasure "$aCGND1CE1";
	set aND1CE1HE1 [measure angle [list $ND1 $CE1 $HE1]];  
	puts   $anglemeasure "$aND1CE1HE1";
	set aND1CE1NE2 [measure angle [list $ND1 $CE1 $NE2]];  
	puts   $anglemeasure "$aND1CE1NE2";
	set aCE1NE2CD2 [measure angle [list $CE1 $NE2 $CD2]];  
	puts   $anglemeasure "$aCE1NE2CD2";
	set aNE2CD2HD2 [measure angle [list $NE2 $CD2 $HD2]];  
	puts   $anglemeasure "$aNE2CD2HD2";

	for {set j 0} {$j <15} {incr j} {
	puts   $anglemeasure "0";
	}
	
	set aCBCAC [measure angle [list $CB $CA $C]];  
	puts   $impranglemeasure "$aCBCAC";
    set dCNCACB [measure dihed [list $C  $N  $CA $CB]]; 
	puts   $dihdmeasure "$dCNCACB";
	set dNCACBHB1 [measure dihed [list $N $CA $CB $HB1]]; 
	puts   $dihdmeasure "$dNCACBHB1";
	set dNCACBHB2 [measure dihed [list $N $CA $CB $HB2]]; 
	puts   $dihdmeasure "$dNCACBHB2";
	set dNCACBCG [measure dihed [list $N $CA $CB $CG]]; 
	puts   $dihdmeasure "$dNCACBCG";
	set dCACBCGND1 [measure dihed [list $CA $CB $CG $ND1]]; 
	puts   $dihdmeasure "$dCACBCGND1";
	set dCBCGND1HD1 [measure dihed [list $CB $CG $ND1 $HD1]]; 
	puts   $dihdmeasure "$dCBCGND1HD1";
	set dCBCGND1CE1 [measure dihed [list $CB $CG $ND1 $CE1]]; 
	puts   $dihdmeasure "$dCBCGND1CE1";
	set dCGND1CE1HE1 [measure dihed [list $CG $ND1 $CE1 $HE1]]; 
	puts   $dihdmeasure "$dCGND1CE1HE1";
	set dCGND1CE1NE2 [measure dihed [list $CG $ND1 $CE1 $NE2]]; 
	puts   $dihdmeasure "$dCGND1CE1NE2";
	set dND1CE1NE2CD2 [measure dihed [list $ND1 $CE1 $NE2 $CD2]]; 
	puts   $dihdmeasure "$dND1CE1NE2CD2";
	set dCE1NE2CD2HD2 [measure dihed [list $CE1 $NE2 $CD2 $HD2]]; 
	puts   $dihdmeasure "$dCE1NE2CD2HD2";
	
	for {set j 0} {$j <15} {incr j} {	
	puts   $dihdmeasure "0";
    }
	}
}

close $bondmeasure
close $impranglemeasure
close $anglemeasure
close $dihdmeasure
close $resi