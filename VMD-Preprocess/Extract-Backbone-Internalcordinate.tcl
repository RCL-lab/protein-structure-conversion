####################################################################################################
## FileName: Extract-Backbone-InternalCoordinate.tcl
## Description: This module extract the internal coordinate of protein backbone
##				Load protein.pdb input in VMD 
##				source Extract-Backbone-InternalCoordinate.tcl
##				Type the number of residue of your protein as input
##				The output are 4 files "BOND.txt","ANGLE.txt","DIHD.txt" 
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
set selatomO [[atomselect top "name O"] get index ]




set bondmeasure [open "BOND.txt" "w"]
for {set i 0} {$i < $NoRes} {incr i} {
set bNCA [measure bond [list [lindex $selatomN $i] [lindex $selatomCA $i]]]; 
puts -nonewline $bondmeasure "NCA"; 
puts -nonewline $bondmeasure " $i:"; 
puts   $bondmeasure "$bNCA";

set bCAC [measure bond [list [lindex $selatomCA $i] [lindex $selatomC $i]]]; 
puts -nonewline $bondmeasure "CAC"; 
puts -nonewline $bondmeasure " $i:"; 
puts   $bondmeasure "$bCAC";

set bCO [measure bond [list [lindex $selatomC $i] [lindex $selatomO $i]]]; 
puts -nonewline $bondmeasure "CO"; 
puts -nonewline $bondmeasure " $i:"; 
puts   $bondmeasure "$bCO"

if {$i< $NoRes - 1} {
set bCN [measure bond [list [lindex $selatomC $i] [lindex $selatomN $i+1]]]; 
puts -nonewline $bondmeasure "CN"; 
puts -nonewline $bondmeasure " $i:"; 
puts   $bondmeasure "$bCN";
};
}; 
close $bondmeasure

set anglemeasure [open "ANGLE.txt" "w"]
for {set i 0} {$i < $NoRes } {incr i} {
set aNCAC [measure angle [list [lindex $selatomN $i] [lindex $selatomCA $i] [lindex $selatomC $i]]]; 
puts -nonewline $anglemeasure "NCAC"; 
puts -nonewline $anglemeasure " $i:"; 
puts   $anglemeasure "$aNCAC";

set aCACO [measure angle [list [lindex $selatomCA $i] [lindex $selatomC $i] [lindex $selatomO $i]]]; 
puts -nonewline $anglemeasure "CACO"; 
puts -nonewline $anglemeasure " $i:"; 
puts   $anglemeasure "$aCACO";

if  { $i< $NoRes - 1} {
set aCACN [measure angle [list [lindex $selatomCA $i] [lindex $selatomC $i] [lindex $selatomN $i+1]]]; 
puts -nonewline $anglemeasure "CACN"; 
puts -nonewline $anglemeasure " $i:"; 
puts   $anglemeasure "$aCACN"

set aCNCA [measure angle [list [lindex $selatomC $i] [lindex $selatomN $i+1] [lindex $selatomCA $i+1]]]; 
puts -nonewline $anglemeasure "CNCA"; 
puts -nonewline $anglemeasure " $i:"; 
puts   $anglemeasure "$aCNCA";
};
}; 
close $anglemeasure

set dihdmeasure [open "DIHD.txt" "w"]
for {set i 0} {$i < $NoRes } {incr i} {
set dNCACO [measure dihed [list [lindex $selatomN $i] [lindex $selatomCA $i] [lindex $selatomC $i] [lindex $selatomO $i]]]; 
puts -nonewline $dihdmeasure "NCACO"; 
puts -nonewline $dihdmeasure " $i:"; 
puts   $dihdmeasure "$dNCACO";

if { $i< $NoRes - 1} {
set dNCACN [measure dihed [list [lindex $selatomN $i] [lindex $selatomCA $i] [lindex $selatomC $i] [lindex $selatomN $i+1]]]; 
puts -nonewline $dihdmeasure "NCACN"; 
puts -nonewline $dihdmeasure " $i:"; 
puts   $dihdmeasure "$dNCACN";

set dCACNCA [measure dihed [list [lindex $selatomCA $i] [lindex $selatomC $i] [lindex $selatomN $i+1] [lindex $selatomCA $i+1]]]; 
puts -nonewline $dihdmeasure "CACNCA"; 
puts -nonewline $dihdmeasure " $i:"; 
puts   $dihdmeasure "$dCACNCA"

set dCNCAC [measure dihed [list [lindex $selatomC $i] [lindex $selatomN $i+1] [lindex $selatomCA $i+1] [lindex $selatomC $i+1]]]; 
puts -nonewline $dihdmeasure "CNCAC"; 
puts -nonewline $dihdmeasure " $i:"; 
puts   $dihdmeasure "$dCNCAC";
};
}; 
close $dihdmeasure
