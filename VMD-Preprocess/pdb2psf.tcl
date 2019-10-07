####################################################################################################
## FileName:    pdb2psf.tcl
## Description: This module converts one PDB file to PSF file.
##				modify the file path of protein.pdb in the module.
##				pdb "/home/VMD/input/5rsanowater.pdb"
##				modify the path where your topology file is located 
##              topology /home/VMD/top_all22_prot_cmap.inp
##				modify the path of psf and pdb files at the end of the module.
##				coordpdb /home/VMD/input/5rsanowater.pdb A 
##              writepdb /home/VMD/output/5rsanowater.pdb
##              writepsf /home/VMD/output/5rsanowater.psf
##				Finally:  source pdb2psf.tcl
##				The output is protein.psf file"
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


resetpsf
package require psfgen

##Modify the path where your topology file is located 
topology /home/VMD/top_all22_prot_cmap.inp
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
segment A {
pdb "/home/VMD/input/5rsanowater.pdb"
first nter
last cter
}
pdbalias atom ILE CD1 CD
pdbalias atom GLY HA2 HA1
pdbalias atom GLY HA3 HA2
pdbalias atom ALA H HN
pdbalias atom ARG H HN
pdbalias atom ARG HB2 HB1
pdbalias atom ARG HB3 HB2
pdbalias atom ARG HG2 HG1
pdbalias atom ARG HG3 HG2
pdbalias atom ARG HD2 HD1
pdbalias atom ARG HD3 HD2
pdbalias atom ASP H HN
pdbalias atom ASP HB2 HB1
pdbalias atom ASP HB3 HB2
pdbalias atom ASN H HN
pdbalias atom ASN HB2 HB1
pdbalias atom ASN HB3 HB2
pdbalias atom CYS HB2 HB1
pdbalias atom CYS HB3 HB2
pdbalias atom CYS H HN
pdbalias atom CYS D HG1
pdbalias atom GLU H HN
pdbalias atom HSD H HN
pdbalias atom HSD HB2 HB1
pdbalias atom HSD HB3 HB2
pdbalias atom HIS H HN
pdbalias atom HIS HB2 HB1
pdbalias atom HIS HB3 HB2
pdbalias atom GLU HB2 HB1
pdbalias atom GLU HB3 HB2
pdbalias atom GLU HG2 HG1
pdbalias atom GLU HG3 HG2
pdbalias atom GLN HB2 HB1
pdbalias atom GLN HB3 HB2
pdbalias atom GLN HG2 HG1
pdbalias atom GLN HG3 HG2
pdbalias atom GLN H HN
pdbalias atom GLY H HN
pdbalias atom ILE H HN
pdbalias atom ILE HG12 HG11
pdbalias atom ILE HG13 HG12
pdbalias atom ILE HD11 HD1
pdbalias atom ILE HD12 HD2
pdbalias atom ILE HD13 HD3
pdbalias atom LEU H HN
pdbalias atom LEU HB2 HB1
pdbalias atom LEU HB3 HB2
pdbalias atom LYS H HN
pdbalias atom LYS HB2 HB1
pdbalias atom LYS HB3 HB2
pdbalias atom LYS HG2 HG1
pdbalias atom LYS HG3 HG2
pdbalias atom LYS HD2 HD1
pdbalias atom LYS HD3 HD2
pdbalias atom LYS HE2 HE1
pdbalias atom LYS HE3 HE2
pdbalias atom MET H HN
pdbalias atom MET HB2 HB1
pdbalias atom MET HB3 HB2
pdbalias atom MET HG2 HG1
pdbalias atom MET HG3 HG2
pdbalias atom PHE H HN
pdbalias atom PHE HB2 HB1
pdbalias atom PHE HB3 HB2
pdbalias atom PRO HB2 HB1
pdbalias atom PRO HB3 HB2
pdbalias atom PRO HD2 HD1
pdbalias atom PRO HD3 HD2
pdbalias atom PRO HG2 HG1
pdbalias atom PRO HG3 HG2
pdbalias atom PRO H HN
pdbalias atom SER H HN
pdbalias atom THR H HN
pdbalias atom TRP H HN
pdbalias atom TRP HB2 HB1
pdbalias atom TRP HB3 HB2
pdbalias atom TYR H HN
pdbalias atom TYR HB2 HB1
pdbalias atom TYR HB3 HB2
pdbalias atom VAL H HN
pdbalias atom ILE CD1 CD


## Modify the following path of coordpdb, psf and pdb files

coordpdb /home/VMD/input/5rsanowater.pdb A 
guesscoord
writepdb /home/VMD/output/5rsanowater.pdb
writepsf /home/VMD/output/5rsanowater.psf

