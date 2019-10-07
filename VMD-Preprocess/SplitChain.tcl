####################################################################################################
## FileName:    SplitChain.tcl
## Description: This module is needed, if the protein that is wanted to be processed consists of multiple chains
##				this module splits the protein to split multiple chains.
##				modify the directory path of pdb file in the module.
##				mol load pdb /home/VMD/input/6awb.pdb
##				modify the path where you want to write the chains pdb
##              $sel writepdb /home/VMD/output/6awb_chain${chain}.pdb
##				Finally:  source SplitChain.tcl
##				The output are different number of chains.
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


##modify the directory path of pdb file in the module.
mol load pdb /home/VMD/input/6awb.pdb
#set water [atomselect top water]
#$water writepdb /home/VMD/output/6awb-water.pdb
set protein [atomselect top protein]
set chains [lsort -unique [$protein get pfrag]]
foreach chain $chains {
set sel [atomselect top "pfrag $chain"]
##	modify the path where you want to write the chains pdb
$sel writepdb /home/VMD/output/6awb_chain${chain}.pdb
}