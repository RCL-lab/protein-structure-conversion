/*####################################################################################################
## FileName:    writeToPDB.cpp
## Description: The module write the output of internal to Cartesian conversion into pdb similar format.
##				The module maps our program output index to the pdb index and omits the dummy atoms from the side chain of each residue.
##				To run this simpily call this function in the either GPU or CPU code version. The two input arguments are backbone Cartesian 
##				coordinate matrix and the side chain Cartesian coordinate matrix. The header file related to this file is ItoC.h.
##				The output is simplepdb.txt.
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
*/

#include <stdlib.h>
#include <stdio.h>                         
#include <iostream>
#include <cstdlib>                                              
#include <time.h>
#include<fstream>
#include <string>
#include <string.h>
#include <iomanip> 
#include <math.h>
#include "ItoC.h"

void protein_file(double* finalxyz,double* sidexyz){
   ofstream mypdb ("simple_pdb.txt");
   ifstream name;
   name.open("resname.txt");
   name.seekg(0, ios::beg);
   name.clear();
   string PRO= "PRO"; char PRO_atom[][5]= {"N", "CA", "C", "O", "CB", "CG","CD", "HA", "HB2", "HB3", "HG2", "HG3", "HD2","HD3","end"};
   string VAL= "VAL"; char VAL_atom[][5]= {"N", "CA", "C", "O", "CB", "CG1", "CG2", "H", "HA", "HB", "HG11", "HG12", "HG13", "HG21","HG22", "HG23","end"};
   string THR= "THR"; char THR_atom[][5]= {"N", "CA", "C", "O", "CB", "OG1", "CG2", "H", "HA", "HB", "HG11", "HG21", "HG22", "HG23","end"};
   string SER= "SER"; char SER_atom[][5]= {"N", "CA", "C", "O", "CB", "OG", "H", "HA", "HB2", "HB3", "HG","end"}; 
   string MET= "MET"; char MET_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "SD", "CE", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HE1", "HE2","HE3","end"};
   string LYS= "LYS"; char LYS_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2","HD3", "HE2", "HE3", "HZ1", "HZ2", "HZ3","end"};
   string LEU= "LEU"; char LEU_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "H", "HA", "HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23","end"};
   string ILE= "ILE"; char ILE_atom[][5]= {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "H", "HA", "HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12","HD13","end"};
   string GLY= "GLY"; char GLY_atom[][5]= {"N", "CA", "C", "O", "H", "HA2", "HA3","end"};
   string GLN= "GLN"; char GLN_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HE21","HE22","end"};
   string GLU= "GLU"; char GLU_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "H", "HA", "HB2", "HB3", "HG2", "HG3","end"};
   string CYS= "CYS"; char CYS_atom[][5]= {"N", "CA", "C", "O", "CB", "CS", "H", "HA", "HB2", "HB3","end"};
   string ASN= "ASN"; char ASN_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "H", "HA", "HB2", "HB3", "HD21", "HD22","end"};
   string ARG= "ARG"; char ARG_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2","HD3", "HE", "HH11", "HH12", "HH13", "HH21", "HH22", "HH23","end"};
   string ALA= "ALA"; char ALA_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD", "HA", "HB2", "HB3", "HG2", "HG3", "HD2","HD3","end"};
   string HIS= "HIS"; char HIS_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1","end"};
   string TRP= "TRP"; char TRP_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2", "H", "HA", "HB2", "HB3", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2","end"};
   string PHE= "PHE"; char PHE_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1","HE2", "HZ","end"};
   string TYR= "TYR"; char TYR_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HH","end"};
   string ASP= "ASP"; char ASP_atom[][5]= {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "H", "HA", "HB2", "HB3","end"};
   

   string str;
   int idxatom=1, idxres=1;
   if (mypdb.is_open() && name.is_open()){
    while(!name.eof()){
        getline(name,str);
		if (str.compare(PRO)==0){
			for (int s=0; s<14; s++){
				if (s<4){
					mypdb<< idxatom<<"  "<<  PRO_atom[s]  <<"  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(PRO_atom[s],"H")==0 || strcmp(PRO_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< PRO_atom[s]  << "  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
				else { 
					if (4<=s&& s<7){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  << "  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<<"\n";
						idxatom++;
				    }
					//when HB,HD,HG ADDED
                    /*if (9<=s && s<=10){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  << PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-8)*3]<<" "<<sidexyz[(idxres-1)*26+(s-8)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-8)*3+2];
	                    idxatom++;
					}
                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  << PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-8)*3+2];
						idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  << PRO <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26]<<" "<<sidexyz[(idxres-1)*26+1]<<" "<< sidexyz[(idxres-1)*26+2];
						idxatom++;
					}*/
					if (9<=s && s<=10){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  << "  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " <<  "na na na" << "\n";
	                    idxatom++;
					}
                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  <<"  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " <<  "na na na" << "\n";
						idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< PRO_atom[s]  <<"  "<< PRO <<"  "<< "A"  <<"  " << idxres  <<"  " <<  "na na na" << "\n";
						idxatom++;
					}
				}
			}
            idxres++;
        }


		if (str.compare(VAL)==0){
            for (int s=0; s<16; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  VAL_atom[s]  <<"  "<< VAL <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(VAL_atom[s],"H")==0 || strcmp(VAL_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<6){
						mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
	                if (s==6){
						mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+s*3]<<" "<<sidexyz[(idxres-1)*26+s*3+1]<<" "<< sidexyz[(idxres-1)*26+s*3+2]<< "\n";
						idxatom++;
					}

                    if (9<=s && s<=12){
						mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-6)*3]<<" "<<sidexyz[(idxres-1)*26+(s-6)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-6)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==15){
						mypdb<< idxatom<<"  "<< VAL_atom[s]  <<"  "<< VAL<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
	            }
		}
        idxres++;
    }

		if (str.compare(THR)==0){
			for (int s=0; s<14; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  THR_atom[s]  <<"  "<< THR <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(THR_atom[s],"H")==0 || strcmp(THR_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<6){
						mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
	                if (s=6){
						mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+s*3]<<" "<<sidexyz[(idxres-1)*26+s*3+1]<<" "<< sidexyz[(idxres-1)*26+s*3+2]<< "\n";
						idxatom++;
					}

                    if (9<=s && s<=10){
						mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==13){
						mypdb<< idxatom<<"  "<< THR_atom[s]  <<"  "<< THR<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-1)*3+2]<< "\n";
						idxatom++;
					}
	            }
			}
        idxres++;
	}

		if (str.compare(SER)==0){
            for (int s=0; s<10; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  SER_atom[s]  <<"  "<< SER <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(SER_atom[s],"H")==0 || strcmp(SER_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< SER_atom[s]  <<"  "<< SER <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                 if (4<=s && s<6){
						mypdb<< idxatom<<"  "<< SER_atom[s]  <<"  "<< SER <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
						idxatom++;
					 }
                     if (8<=s<=9){
						mypdb<< idxatom<<"  "<< SER_atom[s]  <<"  "<< SER <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
	                    idxatom++;
						
					 }
				}
			}
            idxres++;
		}
		
		
		if (str.compare(MET)==0){
            for (int s=0; s<17; s++){
                 if (s<4){
					mypdb<< idxatom<<"  "<<  MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(MET_atom[s],"H")==0 || strcmp(MET_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<7){
						mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+((s-4)*6)*3]<<" "<<sidexyz[(idxres-1)*26+((s-4)*6)*3+1]<<" "<< sidexyz[(idxres-1)*26+((s-4)*6)*3+2]<< "\n";
						idxatom++;
					}
	                if (s==7){
						mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+8)*3]<<" "<<sidexyz[(idxres-1)*26+(s+8)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+8)*3+2]<< "\n";
						idxatom++;
					}

                    if (10<=s && s<=11){
						mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-9)*3]<<" "<<sidexyz[(idxres-1)*26+(s-9)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-9)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (12<=s && s<=13){
						mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-5)*3]<<" "<<sidexyz[(idxres-1)*26+(s-5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-5)*3+2]<< "\n";
						idxatom++;
					}
                    if (14<=s && s<=16){
						mypdb<< idxatom<<"  "<< MET_atom[s]  <<"  "<< MET <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
        idxres++;
		}

		if (str.compare(LYS)==0){
            for (int s=0; s< 22; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(LYS_atom[s],"H")==0 || strcmp(LYS_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<7){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+((s-4)*6)*3]<<" "<<sidexyz[(idxres-1)*26+((s-4)*6)*3+1]<<" "<< sidexyz[(idxres-1)*26+((s-4)*6)*3+2]<< "\n";
						idxatom++;
					}
	                if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3*3+2]<< "\n";
						idxatom++;
					}
                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-6)*3]<<" "<<sidexyz[(idxres-1)*26+(s-6)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-6)*3+2]<< "\n";
						idxatom++;
					}
                    if (15<=s && s<=16){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
						idxatom++;
					}

                    if (17<=s && s<=18){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-1)*3+2]<< "\n";
						idxatom++;
					}
                    if (19<=s && s<=20){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+s*3]<<" "<<sidexyz[(idxres-1)*26+s*3+1]<<" "<< sidexyz[(idxres-1)*26+s*3+2]<< "\n";
						idxatom++;
					}
                    if (s==21){
						mypdb<< idxatom<<"  "<< LYS_atom[s]  <<"  "<< LYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
						idxatom++;
					}
				}
	        }
        idxres++;
	    }

		if (str.compare(LEU)==0){
            for (int s=0; s < 13; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  LEU_atom[s]  <<"  "<< LEU <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(LEU_atom[s],"H")==0 || strcmp(LEU_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s<=7 && s!=6){
						mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
						idxatom++;
					 }
	                if (s==6){
						mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
						idxatom++;
					}

                    if (13<=s && s<=15){
						mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (16<=s && s<=18){
						mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==12){
						mypdb<< idxatom<<"  "<< LEU_atom[s]  <<"  "<< LEU<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-5)*3]<<" "<<sidexyz[(idxres-1)*26+(s-5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-5)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
            idxres++;
		}

		if (str.compare(ILE)==0){
            for (int s=0; s<19; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  ILE_atom[s]  <<"  "<< ILE <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(ILE_atom[s],"H")==0 || strcmp(ILE_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (6<=s && s<=7){	
						mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-5)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-5)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-5)*6*3+2]<< "\n";
						idxatom++;
					}
	                if (4<=s && s<=5){
						mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}

                    if (10<=s && s<=13){
						mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-8)*3]<<" "<<sidexyz[(idxres-1)*26+(s-8)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-8)*3+2]<< "\n";
	                     idxatom++;
					}
                    if (14<=s && s<=15){
						mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
						idxatom++;
					}
                    if (16<=s<=18){
						mypdb<< idxatom<<"  "<< ILE_atom[s]  <<"  "<< ILE<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
				}
	        }
			idxres++;
		}

		if (str.compare(GLY)==0){
            for (int s=0;s<7; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  GLY_atom[s]  <<"  "<< GLY <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				 }
				else if (strcmp(GLY_atom[s],"H")==0 || strcmp(GLY_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< GLY_atom[s]  <<"  "<< GLY <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (s==6){
						mypdb<< idxatom<<"  "<< GLY_atom[s]  <<"  "<< GLY <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-6)*3]<<" "<<sidexyz[(idxres-1)*26+(s-6)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-6)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
		idxres++;
		}

		if (str.compare(GLN)==0){
            for (int s=0; s < 17; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  GLN_atom[s]  <<"  "<< GLN <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(GLN_atom[s],"H")==0 || strcmp(GLN_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<=6){
						mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
						idxatom++;
				    }
	                if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(2*s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(2*s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(2*s-1)*3+2]<< "\n";
						idxatom++;
					}

                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-5)*3]<<" "<<sidexyz[(idxres-1)*26+(s-5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-5)*3+2]<< "\n";
						idxatom++;
					}
                    if (15<=s && s<=16){
						mypdb<< idxatom<<"  "<< GLN_atom[s]  <<"  "<< GLN<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
            idxres++;
		}
		if (str.compare(GLU)==0){
            for (int s=0; s <15; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  GLU_atom[s]  <<"  "<< GLU <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(GLU_atom[s],"H")==0 || strcmp(GLU_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< GLU_atom[s]  <<"  "<< GLU <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<=6){
						mypdb<< idxatom<<"  "<< GLU_atom[s]  <<"  "<< GLU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
						idxatom++;
					}
	                if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< GLU_atom[s]  <<"  "<< GLU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(2*s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(2*s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(2*s-1)*3+2]<< "\n";
						idxatom++;
					}

                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< GLU_atom[s]  <<"  "<< GLU <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< GLU_atom[s]  <<"  "<< GLU<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-5)*3]<<" "<<sidexyz[(idxres-1)*26+(s-5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-5)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
            idxres++;
		}

	if (str.compare(CYS)==0){
        for (int s=0; s<10; s++){
            if (s<4){
				mypdb<< idxatom<<"  "<<  CYS_atom[s]  <<"  "<< CYS <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
				idxatom++;
	        }
	        else if (strcmp(CYS_atom[s],"H")==0 || strcmp(CYS_atom[s],"HA")==0){
				mypdb<< idxatom<<"  "<< CYS_atom[s]  <<"  "<< CYS <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
				idxatom++;
	        }
            else { 
	            if (4<=s && s<=5){
					mypdb<< idxatom<<"  "<< CYS_atom[s]  <<"  "<< CYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
					idxatom++;
				 }

                if (8<=s && s<=9){
					mypdb<< idxatom<<"  "<< CYS_atom[s]  <<"  "<< CYS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
	                idxatom++;
				}
			}
		}
			idxres++;
		}

	if (str.compare(ASN)==0){
        for (int s=0; s<14; s++){
            if (s<4){
				mypdb<< idxatom<<"  "<<  ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
				idxatom++;
	        }
	        else if (strcmp(ASN_atom[s],"H")==0 || strcmp(ASN_atom[s],"HA")==0){
				mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
				idxatom++;
	        }
            else { 
	             if (4<=s && s<=5){
					mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
					idxatom++;
				 }
	             if (s==6){
					mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
					idxatom++;
				 }
	            if (s==7){
					mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+5)*3]<<" "<<sidexyz[(idxres-1)*26+(s+5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+5)*3+2]<< "\n";
					idxatom++;
				 }

                if (10<=s && s<=11){
					mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-9)*3]<<" "<<sidexyz[(idxres-1)*26+(s-9)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-9)*3+2]<< "\n";
	                idxatom++;
				}
                if (12<=s && s<=13){
					mypdb<< idxatom<<"  "<< ASN_atom[s]  <<"  "<< ASN<<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
					idxatom++;
				}
	        }
		}
                       idxres++;
	}

		if (str.compare(ARG)==0){
			for (int s=0; s<24; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(ARG_atom[s],"H")==0 || strcmp(ARG_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
					if (4<=s && s<7){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+((s-4)*6)*3]<<" "<<sidexyz[(idxres-1)*26+((s-4)*6)*3+1]<<" "<< sidexyz[(idxres-1)*26+((s-4)*6)*3+2]<< "\n";
						idxatom++;
					}
	                if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3*3+2]<< "\n";
						idxatom++;
					}
	                if (9<=s && s<=10){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(3*s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(3*s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(3*s-7)*3+2]<< "\n";
						idxatom++;
					}
                    if (11<=s && s<=12){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s<=14){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-12)*3]<<" "<<sidexyz[(idxres-1)*26+(s-12)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-12)*3+2]<< "\n";
						idxatom++;
					}
                    if (15<=s && s<=16){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-8)*3]<<" "<<sidexyz[(idxres-1)*26+(s-8)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-8)*3+2]<< "\n";
						idxatom++;
					}

                    if (17<=s && s<=18){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
                    if (20<=s && s<=21){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
						idxatom++;
					}
                    if (22<=s && s<=23){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
						idxatom++;
					}

                    if (s==19){
						mypdb<< idxatom<<"  "<< ARG_atom[s]  <<"  "<< ARG <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
				}
			}	
            idxres++;
		}

		if (str.compare(ALA)==0){
            for (int s=0; s<10; s++){
                if (s<4){
					mypdb<< idxatom<<"  "<<  ALA_atom[s]  <<"  "<< ALA <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(ALA_atom[s],"H")==0 || strcmp(ALA_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< ALA_atom[s]  <<"  "<< ALA <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (s==4){
						mypdb<< idxatom<<"  "<< ALA_atom[s]  <<"  "<< ALA <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
					
                    if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< ALA_atom[s]  <<"  "<< ALA <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-6)*3]<<" "<<sidexyz[(idxres-1)*26+(s-6)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-6)*3+2]<< "\n";
	                    idxatom++;
					}

                    if (s==9){
						mypdb<< idxatom<<"  "<< ALA_atom[s]  <<"  "<< ALA <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
	                    idxatom++;
					}
				}
			}
            idxres++;
		}
	     
	     
		if (str.compare(HIS)==0){
			for (int s=0; s<16; s++){
				if (s<4){
					mypdb<< idxatom<<"  "<<  HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(HIS_atom[s],"H")==0 || strcmp(HIS_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
				else { 
					if (s==4){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
					if (5<=s && s<=6){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
						idxatom++;
					}
					if (s==7){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
						idxatom++;
					}
					if (8<=s && s<=9){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(2*s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(2*s -10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(2*s -10)*3+2]<< "\n";
						idxatom++;
					}
					if (s==12 || s==13){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-11)*3]<<" "<<sidexyz[(idxres-1)*26+(s-11)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-11)*3+2]<< "\n";
						idxatom++;
					}
					if (s==14 || s==16){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-9)*3]<<" "<<sidexyz[(idxres-1)*26+(s-9)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-9)*3+2]<< "\n";
						idxatom++;
					}
					if (s==15){
						mypdb<< idxatom<<"  "<< HIS_atom[s]  <<"  "<< HIS <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
        idxres++;				
		}


		if (str.compare(TRP)==0){
			for (int s=0; s<24; s++){
				if (s<4){
					mypdb<< idxatom<<"  "<<  TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(TRP_atom[s],"H")==0 || strcmp(TRP_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
				else { 
					if (s==4){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
	                if (5<=s && s<=6){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
						idxatom++;
					}
	                if (7<=s && s<=8){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(-3*s+30)*3]<<" "<<sidexyz[(idxres-1)*26+(-3*s+30)*3+1]<<" "<< sidexyz[(idxres-1)*26+(-3*s+30)*3+2]<< "\n";
						idxatom++;
					 }

                    if (s==10 || s==12){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s)*3]<<" "<<sidexyz[(idxres-1)*26+(s)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (s==9){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-1)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==13){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
						idxatom++;
					}
                    if ( s==11 ){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+5)*3]<<" "<<sidexyz[(idxres-1)*26+(s+5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+5)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==16 || s==17){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-15)*3]<<" "<<sidexyz[(idxres-1)*26+(s-15)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-15)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==18 || s==19){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(2*s-31)*3]<<" "<<sidexyz[(idxres-1)*26+(2*s-31)*3+1]<<" "<< sidexyz[(idxres-1)*26+(2*s-31)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==20 || s==21){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(6*s-109)*3]<<" "<<sidexyz[(idxres-1)*26+(6*s-109)*3+1]<<" "<< sidexyz[(idxres-1)*26+(6*s-109)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==22 || s==23){
						mypdb<< idxatom<<"  "<< TRP_atom[s]  <<"  "<< TRP <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(2*s-31)*3]<<" "<<sidexyz[(idxres-1)*26+(2*s-31)*3+1]<<" "<< sidexyz[(idxres-1)*26+(2*s-31)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
        idxres++;		
		}


		if (str.compare(PHE)==0){
			for (int s=0; s<20; s++){
				if (s<4){
					mypdb<< idxatom<<"  "<<  PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(PHE_atom[s],"H")==0 || strcmp(PHE_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
				else { 
					if (s==4){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
	                if (5<=s && s<=6){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
						idxatom++;
					}
	                if (s==17 || s==19 || s==15){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
						idxatom++;
					 }

                    if (s==10 || s==8){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (13<=s && s<=14){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-12)*3]<<" "<<sidexyz[(idxres-1)*26+(s-12)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-12)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==9){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-1)*3]<<" "<<sidexyz[(idxres-1)*26+(s-1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-1)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==7){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+5)*3]<<" "<<sidexyz[(idxres-1)*26+(s+5)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+5)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==16){
						mypdb<< idxatom<<"  "<< PHE_atom[s]  <<"  "<< PHE <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
        idxres++;
		}

		if (str.compare(TYR)==0){
			for (int s=0; s<21; s++){
				if (s<4){
					mypdb<< idxatom<<"  "<<  TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(TYR_atom[s],"H")==0 || strcmp(TYR_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << "na na na" << "\n";
					idxatom++;
				}
				else { 
					if (s==4){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-4)*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*3+2]<< "\n";
						idxatom++;
					}
	                if (s==7){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+6)*3]<<" "<<sidexyz[(idxres-1)*26+(s+6)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+6)*3+2]<< "\n";
						idxatom++;
					 }
                    if (s==9){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s+2)*3]<<" "<<sidexyz[(idxres-1)*26+(s+2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+2)*3+2]<< "\n";
	                    idxatom++;
					}

                    if (s==10 || s==11 || s==6 || s==8 ||s==5){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-2)*3]<<" "<<sidexyz[(idxres-1)*26+(s-2)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-2)*3+2]<< "\n";
	                    idxatom++;
					}
                    if (14<=s && s<=15){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-13)*3]<<" "<<sidexyz[(idxres-1)*26+(s-13)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-13)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==16 || s==18){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-11)*3]<<" "<<sidexyz[(idxres-1)*26+(s-11)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-11)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==17){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-3)*3]<<" "<<sidexyz[(idxres-1)*26+(s-3)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-3)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==19){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-7)*3]<<" "<<sidexyz[(idxres-1)*26+(s-7)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-7)*3+2]<< "\n";
						idxatom++;
					}
                    if (s==20){
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
						mypdb<< idxatom<<"  "<< TYR_atom[s]  <<"  "<< TYR <<"  "<< "A"  <<"  " << idxres  <<"  " << sidexyz[(idxres-1)*26+(s-10)*3]<<" "<<sidexyz[(idxres-1)*26+(s-10)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-10)*3+2]<< "\n";
						idxatom++;
					}
				}
			}
        idxres++;		
		}
		

		if (str.compare(ASP)==0){
            for (int s=0; s<12; s++){
                 if (s<4){
					mypdb<< idxatom<<"  "<<  ASP_atom[s]  <<"  "<< ASP <<"  "<< "A"  <<"  " << idxres  <<"  " << finalxyz[(idxres-1)*12+s*3] <<" " << finalxyz[(idxres-1)*12+s*3+1]<< " " << finalxyz[(idxres-1)*12+s*3+2]<< "\n";
					idxatom++;
				}
				else if (strcmp(ASP_atom[s],"H")==0 || strcmp(ASP_atom[s],"HA")==0){
					mypdb<< idxatom<<"  "<< ASP_atom[s]  <<"  "<< ASP <<"  "<< "A"  <<"  " << idxres <<"  "  << "na na na" << "\n";
					idxatom++;
				}
                else { 
	                if (4<=s && s<=5){
						mypdb<< idxatom<<"  "<< ASP_atom[s]  <<"  "<< ASP <<"  "<< "A" <<"  " << idxres <<"  " << sidexyz[(idxres-1)*26+(s-4)*6*3]<<" "<<sidexyz[(idxres-1)*26+(s-4)*6*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-4)*6*3+2]<< "\n";
						idxatom++;
					}
	                if (6<=s && s<=7){
						mypdb<< idxatom<<"  "<< ASP_atom[s]  <<"  "<< ASP <<"  "<< "A" <<"  "<< idxres <<"  " << sidexyz[(idxres-1)*26+(s+1)*3]<<" "<<sidexyz[(idxres-1)*26+(s+1)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s+1)*3+2]<< "\n";
						idxatom++;
					}

                    if (10<=s && s<=11){
						mypdb<< idxatom<<"  "<< ASP_atom[s]  <<"  "<< ASP <<"  "<< "A" <<"  "<< idxres <<"  " << sidexyz[(idxres-1)*26+(s-9)*3]<<" "<<sidexyz[(idxres-1)*26+(s-9)*3+1]<<" "<< sidexyz[(idxres-1)*26+(s-9)*3+2]<< "\n";
	                    idxatom++;
					}
                }
			}
            idxres++;    
        }
    }
  }


}
