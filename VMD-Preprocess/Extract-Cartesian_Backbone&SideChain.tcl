####################################################################################################
## FileName:    Extract-Cartesian-Backbone&SideChain.tcl
## Description: This module extract the Cartesian coordinate of both protein Backbone and Sidechain
##				Load protein.pdb input in VMD 
##				source Extract-Cartesian-Backboone&SideChain.tcl
##				The output are 2 files "backbone_xyz.txt", "side_xyz.txt"
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



set sidechain [atomselect top "not backbone"]
set side_xyz [open "side_xyz.txt" "w"]
set sxyz [$sidechain get {x y z}];
foreach i $sxyz j [$sidechain get name] {
 puts  -nonewline $side_xyz "$j"
 puts   -nonewline $side_xyz ":"
 puts   $side_xyz "$i"
 incr   $i
 incr  $j
 }
 close $side_xyz
 
set bb [atomselect top "backbone"]
set bb_xyz [open "backbone_xyz.txt" "w"]
set bxyz [$bb get {x y z}];
foreach i $bxyz j [$bb get name] {
 puts  -nonewline $bb_xyz "$j"
 puts   -nonewline $bb_xyz ":"
 puts   $bb_xyz "$i"
 incr   $i
 incr  $j
 }
close $bb_xyz