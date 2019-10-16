/*####################################################################################################
## FileName:    InternaltoCartesian_CPU.cpp
## Description: The module read the intenal coordinate of a protein and serially convert it to cartesian coordinate by walking on the protein chain.
##				This is the serial CPU version of the program. The input files of this program is Bond, Angle, Dihd, imprANGLE_SIDE,BOND_SIDE, ANGLE_SIDE
##				DIHD_SIDE,resname. the output can be two raw data file, backboneXYZ, sideXYZ. 
##				by including ItoC.h header file, we can call writeToPDB.cpp and get the pdb like output. (you only need to uncomment that part)
##				The output is simplepdb.txt.
##			    To run this file just need to compile it by g++ and then run it:  
## 				g++ InternalToCartesian_CPU.cpp -o InternalToCatesian_out
##				./InternalToCatesian_out
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
#include <fstream>
#include <string>
#include <iomanip> 
#include <math.h>
#include "ItoC.h"
#include <sys/time.h>

using namespace std; 

void findXYZ_Dihedral(double *A, double *B, double *C, double CD, double BCD, double ABCD,double *D);
double *temproryXYZ(struct bond *batom,struct angle *angleaatom,struct dihd *dihatom,int noresidue, double *N_P, double *CA_P, double *C_P);
double *finalXYZ(struct bond *batom,struct angle *angleaatom,struct dihd *dihatom, double *tempXYZ,int noresidue);
double  norm(struct XYZ A);
struct XYZ cross(struct XYZ A,struct XYZ B);
void findXYZ_ImproperDihedralTemp(double *A, double *B, double *D, double AC,double BAC, double ABCD,double *C);
void findXYZ_ImproperDihedral(double *A, double *B, double *D, double AC, double DAC,double BAC, double ABCD,double *C);
double *computeSideXYZ(double *bside, double *angleside, double *dihside, double *angleimprside, int noresidue, double *globalXYZ);
double getTorsion (double *A, double *B, double *C, double *D);


double *calcrotationmatrix(double BCD,double *N);
double *mergrot(double *v1,double *v2,double *v3);
double  norm(double *A);
double* Trans(double *A,int s);
double *reg_matrixmult(double *A, double *B, int m, int n,int p);
double *reg_cross(double *A,double *B);
double *buildup( double *tempXYZ,int noresidue, double* Translation);
double *matrixmult(double *thetarotation,double *D);



#define pi          3.141592653589793238462643383279502884L


/*struct bond{
	double NCA;
	double CAC;
	double CO;
	double CN;
};

struct angle{
	double NCAC;
	double CACO;
	double CACN;
	double CNCA;
};

struct dihd{
	double NCACO;
	double NCACN;
	double CACNCA;
	double CNCAC;
};*/

struct XYZ{
	double X;
	double Y;
	double Z;
};

int main()
{
	
  ifstream mybond,myangle,mydih, mybond_side, myangle_side,myangle_imprside, mydih_side;
  mybond.open ("B_5rsa_noindex.txt");
  myangle.open ("A_5rsa_noindex.txt");
  mydih.open ("DIHD_5rsa_noindex.txt");

  mybond_side.open ("BOND_SIDE.txt");
  myangle_side.open ("ANGLE_SIDE.txt");
  myangle_imprside.open ("imprANGLE_SIDE.txt");
  mydih_side.open ("DIHD_SIDE.txt");


  string  str;
  int line=0;
  int i=0;
  int atomperresidue=26;
    
  while (!mybond.eof()){
	getline(mybond,str); 
	line++;
   }

   int noresidue=line/4;
   bond *batom=(bond*) malloc(noresidue*sizeof(bond));
   angle *angleatom=(angle*)malloc(noresidue*sizeof(angle));
   dihd *dihatom=(dihd*)malloc(noresidue*sizeof(dihd));
   double *tempXYZ=(double*)malloc(noresidue*12*sizeof(double));

   double *bside=(double*) malloc(noresidue*atomperresidue*sizeof(double));
   double *angleside=(double*)malloc(noresidue*atomperresidue*sizeof(double));
   double *angleimprside=(double*)malloc(noresidue*sizeof(double));
   double *dihside=(double*)malloc(noresidue*atomperresidue*sizeof(double));
   double *out_sideXYZ=(double*)malloc(noresidue*atomperresidue*3*sizeof(double));



/////////////////////////////////////////////////
   mybond.clear();
   mybond.seekg(0, ios::beg);

   while (!mybond.eof()){
		getline(mybond,str); 
		batom[i].NCA=atof(str.c_str());

		getline(mybond,str); 
		batom[i].CAC=atof(str.c_str());

		getline(mybond,str); 
		batom[i].CO=atof(str.c_str());

		getline(mybond,str); 
		batom[i].CN=atof(str.c_str());

		i++;

	}
  mybond.close();

		i=0;
	    while (!myangle.eof()){

		getline(myangle,str); 
		angleatom[i].NCAC=atof(str.c_str());

		getline(myangle,str); 
		angleatom[i].CACO=atof(str.c_str());

		getline(myangle,str); 
		angleatom[i].CACN=atof(str.c_str());

		getline(myangle,str); 
		angleatom[i].CNCA=atof(str.c_str());

		i++;

	}


  myangle.close();		
        i=0;
	    while (!mydih.eof()){
		getline(mydih,str); 
		dihatom[i].NCACO=atof(str.c_str());

		getline(mydih,str); 
		dihatom[i].NCACN=atof(str.c_str());

		getline(mydih,str); 
		dihatom[i].CACNCA=atof(str.c_str());

		getline(mydih,str); 
		dihatom[i].CNCAC=atof(str.c_str());

		i++;

	}
    mydih.close();	

   // reading the sidechain input
    int i1=0;
    while (i1<noresidue*26){
		getline(mybond_side,str); 
		bside[i1]=atof(str.c_str());

		getline(myangle_side,str); 
		angleside[i1]=atof(str.c_str());

		getline(mydih_side,str); 
		dihside[i1]=atof(str.c_str());
		i1++;
	}
     i1=0;
     while (!myangle_imprside.eof()){
              getline(myangle_imprside,str); 
              angleimprside[i1]=atof(str.c_str());
              i1++;
     }


  mybond_side.close();
  myangle_side.close();
  mydih_side.close();

  //////////////////////////////////////////////////////
 
  double N_P[3]={18.39299964904785,24.60700035095215,0.7289999723434448};
  double CA_P[3]={19.19499969482422,24.738000869750977,1.9490000009536743};
  double C_P[3]={20.62299919128418,25.194000244140625,1.684000015258789};


  tempXYZ   = temproryXYZ(batom,angleatom,dihatom,noresidue, N_P, CA_P,C_P);
   
  free (batom);
  free (angleatom);
  free (dihatom);
  
  out_sideXYZ   = computeSideXYZ(bside, angleside, dihside, angleimprside, noresidue,  tempXYZ); 
  ofstream myfile ("backboneXYZ.txt");
   cout<< noresidue<<"\n";
   if (myfile.is_open())
   {

	   	   for(int k=0;k<noresidue;k++){
			  
			   myfile<<"N: "<<tempXYZ[k*12]<<","<<tempXYZ[k*12+1]<<","<<tempXYZ[k*12+2]<<"\n";
			   myfile<<"CA: "<<tempXYZ[k*12+3]<<","<<tempXYZ[k*12+4]<<","<<tempXYZ[k*12+5]<<"\n";
			   myfile<<"C: "<<tempXYZ[k*12+6]<<","<<tempXYZ[k*12+7]<<","<<tempXYZ[k*12+8]<<"\n";
			   myfile<<"O: "<<tempXYZ[k*12+9]<<","<<tempXYZ[k*12+10]<<","<<tempXYZ[k*12+11]<<"\n";

	   }

    myfile.close();
   }
   
  
   else cout << "Unable to open file";


   ofstream myside ("sideXYZ.txt");

   if (myside.is_open())
   {
	for(int k=0;k<noresidue;k++){
		for (int a = 0; a<26 ; a++) {
			myside<< out_sideXYZ[k*26+a*3] << " " <<out_sideXYZ[k*26+a*3+1] << " " <<out_sideXYZ[k*26+a*3+2]<<"\n";

		   }
        }
			   
    

    myside.close();
   }
  else cout << "Unable to open file";
  
  //writeTopdb file output, you can uncomment it to get a combined pdb like output// 
  /* protein_file(tempXYZ,out_sideXYZ); */
  
  return 0;

}



void findXYZ_Dihedral(double *A, double *B, double *C, double CD, double BCD, double ABCD,double *D){

	double IJ[3],JK[3],N[3],DC[3],BC[3],CrossResult[3],rotationmatrix[9];
       for(int i=0;i<3;i++){
		IJ[i]=B[i]-A[i];
		BC[i]=C[i]-B[i];
       }
    
	double normJK=sqrt(pow(BC[0],2)+pow(BC[1],2)+pow(BC[2],2));

       for (int i=0;i<3;i++){
		JK[i]=BC[i]/normJK;
       }
       
       /*CROSS IJ,JK*/
	CrossResult[0]=IJ[1]*JK[2]-IJ[2]*JK[1];
	CrossResult[1]=-(IJ[0]*JK[2]-IJ[2]*JK[0]);
	CrossResult[2]=IJ[0]*JK[1]-IJ[1]*JK[0];
       /*End CROSS IJ,JK*/ 
       double normCross=sqrt(pow(CrossResult[0],2)+pow(CrossResult[1],2)+pow(CrossResult[2],2));
       
	for(int i=0;i<3;i++){
	N[i]=CrossResult[i]/normCross;
	D[i]=C[i]+(CD*JK[i]);
                 DC[i]=D[i]-C[i];
       }

       double pBCD=(pi-(pi*BCD/180));
       
	/*Calrotationmatrix*/
	double s=sinf(pBCD);
	double c=cosf(pBCD);
	rotationmatrix[0]=pow(N[0],2)*(1-c)+c;
	rotationmatrix[1]=N[0]*N[1]*(1-c)-N[2]*s;
	rotationmatrix[2]=N[0]*N[2]*(1-c)+N[1]*s;
	rotationmatrix[3]=N[0]*N[1]*(1-c)+N[2]*s;
	rotationmatrix[4]=pow(N[1],2)*(1-c)+c;
	rotationmatrix[5]=N[1]*N[2]*(1-c)-N[0]*s;
	rotationmatrix[6]=N[0]*N[2]*(1-c)-N[1]*s;
	rotationmatrix[7]=N[2]*N[1]*(1-c)+N[0]*s;
	rotationmatrix[8]=pow(N[2],2)*(1-c)+c;
	/*End Calrotationmatrix*/
	
	/*Matrix multiplication*/	  
	D[0]=rotationmatrix[0]*DC[0]+rotationmatrix[1]*DC[1]+rotationmatrix[2]*DC[2];
	D[1]=rotationmatrix[3]*DC[0]+rotationmatrix[4]*DC[1]+rotationmatrix[5]*DC[2];
	D[2]=rotationmatrix[6]*DC[0]+rotationmatrix[7]*DC[1]+rotationmatrix[8]*DC[2];
       /*End Matrixmultiplication*/   


	for(int i=0;i<3;i++){ 
       D[i]=D[i]+C[i];
	}

	double pABCD=((pi*ABCD)/180);
	
	/*Calrotationmatrix*/
	s=sinf(pABCD);
	c=cosf(pABCD);
	rotationmatrix[0]=pow(JK[0],2)*(1-c)+c;
	rotationmatrix[1]=JK[0]*JK[1]*(1-c)-JK[2]*s;
	rotationmatrix[2]=JK[0]*JK[2]*(1-c)+JK[1]*s;
	rotationmatrix[3]=JK[0]*JK[1]*(1-c)+JK[2]*s;
	rotationmatrix[4]=pow(JK[1],2)*(1-c)+c;
	rotationmatrix[5]=JK[1]*JK[2]*(1-c)-JK[0]*s;
	rotationmatrix[6]=JK[0]*JK[2]*(1-c)-JK[1]*s;
	rotationmatrix[7]=JK[2]*JK[1]*(1-c)+JK[0]*s;
	rotationmatrix[8]=pow(JK[2],2)*(1-c)+c;
	/*End Calrotationmatrix*/
       
	for(int i=0;i<3;i++)
       	DC[i]=D[i]-C[i];
       

	/*Matrix multiplication*/	  
	D[0]=rotationmatrix[0]*DC[0]+rotationmatrix[1]*DC[1]+rotationmatrix[2]*DC[2];
	D[1]=rotationmatrix[3]*DC[0]+rotationmatrix[4]*DC[1]+rotationmatrix[5]*DC[2];
	D[2]=rotationmatrix[6]*DC[0]+rotationmatrix[7]*DC[1]+rotationmatrix[8]*DC[2];
       /*End Matrixmultiplication*/  

       for(int i=0;i<3;i++)	
		D[i]=D[i]+C[i];
       
       

}


void findXYZ_ImproperDihedralTemp(double *A, double *B, double *D, double AC,double BAC, double ABCD,double *C)
{
    double IJ[3],JK[3],N[3],CA[3],AD[3],CrossResult[3],rotationmatrix[9];
    for(int i=0;i<3;i++){
		IJ[i]=B[i]-A[i];
		AD[i]=A[i]-D[i];
    }
    
    double normJK=sqrt(pow(AD[0],2)+pow(AD[1],2)+pow(AD[2],2));

    for (int i=0;i<3;i++){
		JK[i]=AD[i]/normJK;
    }
       
    /*CROSS IJ,JK*/
    CrossResult[0]=IJ[1]*JK[2]-IJ[2]*JK[1];
    CrossResult[1]=-(IJ[0]*JK[2]-IJ[2]*JK[0]);
    CrossResult[2]=IJ[0]*JK[1]-IJ[1]*JK[0];
    /*End CROSS IJ,JK*/ 
    double normCross=sqrt(pow(CrossResult[0],2)+pow(CrossResult[1],2)+pow(CrossResult[2],2));
    for(int i=0;i<3;i++){
		N[i]=CrossResult[i]/normCross;
		C[i]=A[i]+(AC*JK[i]);
       	CA[i]=C[i]-A[i];
     }

       double pBAC=(pi-(pi*BAC/180));
       
       /*Calrotationmatrix*/
       double s=sinf(pBAC);
       double c=cosf(pBAC);
       rotationmatrix[0]=pow(N[0],2)*(1-c)+c;
       rotationmatrix[1]=N[0]*N[1]*(1-c)-N[2]*s;
       rotationmatrix[2]=N[0]*N[2]*(1-c)+N[1]*s;
       rotationmatrix[3]=N[0]*N[1]*(1-c)+N[2]*s;
       rotationmatrix[4]=pow(N[1],2)*(1-c)+c;
       rotationmatrix[5]=N[1]*N[2]*(1-c)-N[0]*s;
       rotationmatrix[6]=N[0]*N[2]*(1-c)-N[1]*s;
       rotationmatrix[7]=N[2]*N[1]*(1-c)+N[0]*s;
       rotationmatrix[8]=pow(N[2],2)*(1-c)+c;
       /*End Calrotationmatrix*/
	
       /*Matrix multiplication*/	  
       C[0]=rotationmatrix[0]*CA[0]+rotationmatrix[1]*CA[1]+rotationmatrix[2]*CA[2];
       C[1]=rotationmatrix[3]*CA[0]+rotationmatrix[4]*CA[1]+rotationmatrix[5]*CA[2];
       C[2]=rotationmatrix[6]*CA[0]+rotationmatrix[7]*CA[1]+rotationmatrix[8]*CA[2];
       /*End Matrixmultiplication*/   


       for(int i=0;i<3;i++){ 
       	C[i]=C[i]+A[i];
        }

        double pABCD=((pi*ABCD)/180);
	
        /*Calrotationmatrix*/
        s=sinf(pABCD);
        c=cosf(pABCD);
        rotationmatrix[0]=pow(JK[0],2)*(1-c)+c;
        rotationmatrix[1]=JK[0]*JK[1]*(1-c)-JK[2]*s;
        rotationmatrix[2]=JK[0]*JK[2]*(1-c)+JK[1]*s;
        rotationmatrix[3]=JK[0]*JK[1]*(1-c)+JK[2]*s;
        rotationmatrix[4]=pow(JK[1],2)*(1-c)+c;
        rotationmatrix[5]=JK[1]*JK[2]*(1-c)-JK[0]*s;
        rotationmatrix[6]=JK[0]*JK[2]*(1-c)-JK[1]*s;
        rotationmatrix[7]=JK[2]*JK[1]*(1-c)+JK[0]*s;
        rotationmatrix[8]=pow(JK[2],2)*(1-c)+c;
        /*End Calrotationmatrix*/
       
       for(int i=0;i<3;i++){
       	CA[i]=C[i]-A[i];
       }

       /*Matrix multiplication*/	  
        C[0]=rotationmatrix[0]*CA[0]+rotationmatrix[1]*CA[1]+rotationmatrix[2]*CA[2];
        C[1]=rotationmatrix[3]*CA[0]+rotationmatrix[4]*CA[1]+rotationmatrix[5]*CA[2];
        C[2]=rotationmatrix[6]*CA[0]+rotationmatrix[7]*CA[1]+rotationmatrix[8]*CA[2];
       /*End Matrixmultiplication*/  

       for(int i=0;i<3;i++){	
		C[i]=C[i]+A[i];
       }
       

}

 void findXYZ_ImproperDihedral(double *A, double *B, double *D, double AC, double DAC,double BAC, double ABCD,double *C){
       

       double n1[3],n2[3],r0[3],n3[3],CrossResult[3],r0A[3];
       double h1n1n2=0,nra=0,T=0,diffABCD=1000,h2n1n2=0,n1n2=0;
       double phi,C1[3],C2[3],T1,T2;
       double phi_start=-0.0005;
       double phi_interval=0.00001;
       double phi_end=0.0005;
	
       ABCD=pi*ABCD/180;
       for(int i=0;i<3;i++){
		n1[i]=B[i]-A[i];
		n2[i]=D[i]-A[i];
       }
    
       double normn1=sqrt(pow(n1[0],2)+pow(n1[1],2)+pow(n1[2],2));
       double normn2=sqrt(pow(n2[0],2)+pow(n2[1],2)+pow(n2[2],2));

       for (int i=0;i<3;i++){
       	n1[i]=n1[i]/normn1;
       	n2[i]=n2[i]/normn2;
       	n1n2+=n1[i]*n2[i];
       }
       
       double h1=AC*cosf(pi*BAC/180)+n1[0]*A[0]+n1[1]*A[1]+n1[2]*A[2];
       double h2=AC*cosf(pi*DAC/180)+n2[0]*A[0]+n2[1]*A[1]+n2[2]*A[2];
       
       double phrase=1-pow(n1n2,2);
       for (int i=0;i<3;i++){
		h2n1n2+=h2*n1[i]*n2[i];
       	h1n1n2+=h1*n1[i]*n2[i];
       }

       double c1=(h1-(h2n1n2))/phrase;
       double c2=(h2-(h1n1n2))/phrase;
       

       for (int i=0;i<3;i++){
       	r0[i]=c1*n1[i]+c2*n2[i];
       }

       /*CROSS n1,n2*/
       CrossResult[0]=n1[1]*n2[2]-n1[2]*n2[1];
       CrossResult[1]=-(n1[0]*n2[2]-n1[2]*n2[0]);
       CrossResult[2]=n1[0]*n2[1]-n1[1]*n2[0];
       /*End CROSS n1,n2*/ 
       double normCross=sqrt(pow(CrossResult[0],2)+pow(CrossResult[1],2)+pow(CrossResult[2],2));
       
       for(int i=0;i<3;i++){
		n3[i]=CrossResult[i]/normCross;   
       	r0A[i]=r0[i]-A[i];
       	nra+=n3[i]*(r0A[i]);
       }
       
       double normr0A=sqrt(pow(r0A[0],2)+pow(r0A[1],2)+pow(r0A[2],2));
       double r_d1,i_d1;

       if((pow(nra,2)-pow(normr0A,2)+pow(AC,2))<0)
       {  
          	i_d1=sqrt(-(pow(nra,2)-pow(normr0A,2)+pow(AC,2)));
          	r_d1=-nra;
        }
       else
       {
          	i_d1=0;
          	r_d1=-nra+sqrt(pow(nra,2)-pow(normr0A,2)+pow(AC,2));
       }
       
       double d2=-nra-sqrt(pow(nra,2)-pow(normr0A,2)+pow(AC,2));
       if(i_d1!=0)
       {
        for(phi=phi_start;phi<=phi_end;phi+=phi_interval)
          {
                findXYZ_ImproperDihedralTemp(A,B,D,AC,BAC,phi,C1);
                T=getTorsion(A,B,C,D);
                if(fabs(T-ABCD) < diffABCD)
                {
                     diffABCD=fabs(T-ABCD);
                     for (int i=0;i<3;i++)
                          C[i]=C1[i];

                  }		
             }              
        }
      else
       {
        for (int i=0;i<3;i++)
        {
             C1[i]=r0[i]+r_d1*n3[i];
             C2[i]=r0[i]+d2*n3[i];
        } 
        T1=getTorsion(A,B,C1,D);
        T2=getTorsion(A,B,C2,D);
        
        if(fabs(T1-ABCD)<fabs(T2-ABCD)){
          for (int i=0;i<3;i++)
              C[i]=C1[i];
       }
        else{
          for (int i=0;i<3;i++)
               C[i]=T1;
          }
        }

}

double getTorsion (double *A, double *B, double *C, double *D)
{
 	double IJ[3],JK[3],KL[3],term1[3],term2[3],term1_1,term2_2,term3[3],T;
       double normJK=0;
 	for(int i=0;i<3;i++)
       {
           IJ[i]=B[i]-A[i];
           JK[i]=C[i]-B[i];
           KL[i]=D[i]-C[i];
           normJK+=pow(JK[i],2);
        }   
    normJK=sqrt(normJK);
  	for(int i=0;i<3;i++)
            term1[i]=normJK*IJ[i];
        
       /*CROSS JK,KL*/
	term2[0]=JK[1]*KL[2]-JK[2]*KL[1];
	term2[1]=-(JK[0]*KL[2]-JK[2]*KL[0]);
	term2[2]=JK[0]*KL[1]-JK[1]*KL[0];
       /*End CROSS JK,KL*/ 

       /*CROSS IJ,JK*/
	term3[0]=IJ[1]*JK[2]-IJ[2]*JK[1];
	term3[1]=-(IJ[0]*JK[2]-IJ[2]*JK[0]);
	term3[2]=IJ[0]*JK[1]-IJ[1]*JK[0];
       /*End CROSS IJ,JK*/ 

    term1_1=term1[0]*term2[0]+term1[1]*term2[1]+term1[2]*term2[2];
    term2_2=term2[0]*term3[0]+term2[1]*term3[1]+term2[2]*term3[2];
       
       
	T=atan2(term1_1,term2_2);
     
       return T;
 
}


double *computeSideXYZ(double *bside, double *angleside, double *dihside, double *angleimprside, int noresidue, double *globalXYZ){


    double bond[26], angle[26], dihd[26];
	double A[3], B[3], C[3], D[3], temp[3];    	
    double AC, BAC, DAC, ABCD, CD, BCD;
	int atomperresidue =26;
    double *sideXYZ= (double*)malloc(noresidue*atomperresidue*3*sizeof(double));
	int tid =0;

  	ifstream resname;
	string PRO= "PRO";
   	string VAL= "VAL";
   	string THR= "THR";
   	string SER= "SER";
   	string MET= "MET";
   	string LYS= "LYS";
   	string LUE= "LUE";
   	string ILE= "ILE";
   	string GLY= "GLY";
   	string GLN= "GLN";
   	string GLU= "GLU";
   	string CYS= "CYS";
   	string ASN= "ASN";
   	string ARG= "ARG";
   	string ALA= "ALA";
   	string HIS= "HIS";
   	string TRP= "TRP";
   	string PHE= "PHE";
   	string TYR= "TYR";
   	string ASP= "ASP";

  	resname.open ("resname.txt");
	string str;
	int ii = 0; 
 	while (!resname.eof()){
 		getline(resname,str);

		if (str.compare(PRO)==0){
				cout<<PRO<<"\n";
			    // compute CB
			    for (int kk=0;kk<3;kk++){
					A[kk]=globalXYZ[(ii+1)*12+3+kk];
					B[kk]=globalXYZ[(ii+1)*12+kk];
					D[kk]=globalXYZ[(ii+1)*12+6+kk];
          		}
			    AC=bside[ii*atomperresidue];
          		BAC=angleside[ii*atomperresidue];
          		DAC=angleimprside[ii];
          		ABCD=dihside[ii*atomperresidue];
	   		    findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
          		for (int kk=0;kk<3;kk++){
          			 sideXYZ[ii*atomperresidue*3+kk]=C[kk];
				     cout<<A[kk]<<", ";
					 cout<<B[kk]<<", ";
					 cout<<D[kk]<<", ";
					 cout<<sideXYZ[ii*atomperresidue*3+kk]<<", ";
	  			}  
				cout <<"\n";
			    // compute CG 
         		for (int kk=0;kk<3;kk++){
          			B[kk]=globalXYZ[(ii+1)*12+3+kk];
          			A[kk]=globalXYZ[(ii+1)*12+kk];   
	  			}
	   			CD=bside[ii*atomperresidue+1];
          		BCD=angleside[ii*atomperresidue+1];
          		ABCD=dihside[ii*atomperresidue+1];          
	   			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);

			    //Compute CD         
          	    for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
	   			} 
          		CD=bside[ii*atomperresidue+2];
                BCD=angleside[ii*atomperresidue+2];
                ABCD=dihside[ii*atomperresidue+2]; 
	            findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		    for (int kk=0; kk<3;  kk++){
					sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
	   			} 

			                //Till we have all the H in VMD processing
          			               /*CD=bside[ii*atomperresidue+3];
          			               BCD=angleside[ii*atomperresidue+3];
          			               ABCD=dihside[ii*atomperresidue+3]; 
          			               findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          			               for (int kk=0;kk<3;kk++){
          				sideXYZ[ii*atomperresidue*3+kk+9]=D[kk];   
          			                } 
         			               for (int kk=0;kk<3;kk++){
          				A[kk]=B[kk];
          				B[kk]=C[kk];
          				C[kk]=D[kk];   
          				}

          			               CD=bside[ii*atomperresidue+4];
                                                                  BCD=angleside[ii*atomperresidue+4];
                                                                  ABCD=dihside[ii*atomperresidue+4];          
	                                                 findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
                                                                  for (int kk=0;kk<3;kk++){
                                                                       sideXYZ[ii*atomperresidue*3+kk+12]=D[kk];   
	                                                  } 
          			               CD=bside[ii*atomperresidue+5];
          			               BCD=angleside[ii*atomperresidue+5];
          			               ABCD=dihside[ii*atomperresidue+5]; 
	   		               findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          			               for (int kk=0;kk<3;kk++){
          				sideXYZ[ii*atomperresidue*3+kk+15]=D[kk];   
	   			} 
	   		               CD=bside[ii*atomperresidue+6];
          			               BCD=angleside[ii*atomperresidue+6];
          			               ABCD=dihside[ii*atomperresidue+6]; 
          			               findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		               for (int kk=0;kk<3;kk++){
          				sideXYZ[ii*atomperresidue*3+kk+18]=D[kk];   
	   			} 
        			              for (int kk=0;kk<3;kk++){
          				A[kk]=B[kk];
          				B[kk]=C[kk];
         				C[kk]=D[kk];   
	  		               }
          			             CD=bside[ii*atomperresidue+7];
          			             BCD=angleside[ii*atomperresidue+7];
          			             ABCD=dihside[ii*atomperresidue+7];          
	   		             findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          			             for (int kk=0;kk<3;kk++){
          				sideXYZ[ii*atomperresidue*3+kk+21]=D[kk];   
	  		              } 
          			             CD=bside[ii*atomperresidue+8];
          			             BCD=angleside[ii*atomperresidue+8];
                                                                ABCD=dihside[ii*atomperresidue+8]; 
                                                                findXYZ_Dihedral(A, B, C, CD,BCD,ABCD,D);
                                                                for (int kk=0;kk<3;kk++){
                                                                     sideXYZ[ii*atomperresidue*3+kk+24]=D[kk];   
	  		             } */
		ii++;
		}
		
        else if (str.compare(HIS)==0){
			    cout << HIS<< "\n";
    			for (int kk=0;kk<0;kk++){
					A[kk]=globalXYZ[(ii+1)*12+3+kk];
					B[kk]=globalXYZ[(ii+1)*12+kk];
					D[kk]=globalXYZ[(ii+1)*12+6+kk];
          		}
          		AC=bside[ii*atomperresidue];
          		BAC=angleside[ii*atomperresidue];
          		DAC=angleimprside[ii];
          		ABCD=dihside[ii*atomperresidue];
				//Calculation CB
          		findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk]=C[kk];   
          		} 
				
				
				//rearrange for calculation HB1,HB2
         		for (int kk=0;kk<3;kk++){
          			B[kk]=globalXYZ[(ii+1)*12+3+kk];
          			A[kk]=globalXYZ[(ii+1)*12+kk];   
          		}
          		CD=bside[ii*atomperresidue+1];
          		BCD=angleside[ii*atomperresidue+1];
          		ABCD=dihside[ii*atomperresidue+1]; 
				//computation HB1
          		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
				} 
          		CD=bside[ii*atomperresidue+2];
          		BCD=angleside[ii*atomperresidue+2];
          		ABCD=dihside[ii*atomperresidue+2];
				//computation HB2
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
				for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
				} 
          		CD=bside[ii*atomperresidue+3];
          		BCD=angleside[ii*atomperresidue+3];
          		ABCD=dihside[ii*atomperresidue+3]; 
				//computation CG
          		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
         		sideXYZ[ii*atomperresidue*3+kk+9]=D[kk];   
				}
				
				
				//rearrange for calculation ND1				
         		for (int kk=0;kk<3;kk++){
          		     A[kk]=B[kk];
          			 B[kk]=C[kk];
          			 C[kk]=D[kk];   
				}

          		CD=bside[ii*atomperresidue+4];
          		BCD=angleside[ii*atomperresidue+4];
          		ABCD=dihside[ii*atomperresidue+4];  
				//computation ND1				
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+12]=D[kk];   
			   } 
			   
			   
				//rearrange for calculation HD1, CE1
         		for (int kk=0;kk<3;kk++){
          			A[kk]=B[kk];
          			B[kk]=C[kk];
          			C[kk]=D[kk];   
				}
          		CD=bside[ii*atomperresidue+5];
          		BCD=angleside[ii*atomperresidue+5];
          		ABCD=dihside[ii*atomperresidue+5]; 
				//computation HD1
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+15]=D[kk];   
				} 
				CD=bside[ii*atomperresidue+6];
          		BCD=angleside[ii*atomperresidue+6];
          		ABCD=dihside[ii*atomperresidue+6]; 
				//computation CE1
          		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
				for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+18]=D[kk];   
				} 
				
				
				//rearrange for calculation HE1,NE2
        		for (int kk=0;kk<3;kk++){
          			A[kk]=B[kk];
          			B[kk]=C[kk];
          			C[kk]=D[kk];   
				}
				CD=bside[ii*atomperresidue+7];
          		BCD=angleside[ii*atomperresidue+7];
          		ABCD=dihside[ii*atomperresidue+7];   
				
				//computation HE1
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+21]=D[kk];   
				} 
          		CD=bside[ii*atomperresidue+8];
          		BCD=angleside[ii*atomperresidue+8];
          		ABCD=dihside[ii*atomperresidue+8]; 
				
				//computation NE2
				findXYZ_Dihedral(A, B, C, CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+24]=D[kk];   
				}
				//rearrange for calculation CD2
        		for (int kk=0;kk<3;kk++){
          				A[kk]=B[kk];
          				B[kk]=C[kk];
          				C[kk]=D[kk];   
				}
				CD=bside[ii*atomperresidue+9];
          		BCD=angleside[ii*atomperresidue+9];
          		ABCD=dihside[ii*atomperresidue+9]; 
				//computation CD2
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          		for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
				}
				
				//rearrange for calculation HD2
        		for (int kk=0;kk<3;kk++){
          			A[kk]=B[kk];
          			B[kk]=C[kk];
          			C[kk]=D[kk];   
				}
          		CD=bside[ii*atomperresidue+10];
          		BCD=angleside[ii*atomperresidue+10];
          		ABCD=dihside[ii*atomperresidue+10];   
				//computation HD2
				findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	    for (int kk=0;kk<3;kk++){
          			sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
				} 
			ii++;
       		}
		

		else if (str.compare(TRP)==0){
			 cout << TRP<< "\n";
        	for (int kk=0;kk<3;kk++){
          		A[kk]=globalXYZ[(ii+1)*12+3+kk];
          		B[kk]=globalXYZ[(ii+1)*12+kk];
          		D[kk]=globalXYZ[(ii+1)*12+6+kk];
	  		}
	   		AC=bside[ii*atomperresidue];
          	BAC=angleside[ii*atomperresidue];
          	DAC=angleimprside[ii];
          	ABCD=dihside[ii*atomperresidue];
			//computing CB
	   		findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk]=C[kk];   
	  		}   
			
			//rearrange for computing CG, HB1,HB2
         	for (int kk=0;kk<3;kk++){
          		B[kk]=globalXYZ[(ii+1)*12+3+kk];
          		A[kk]=globalXYZ[(ii+1)*12+kk];   
	  		}
	   		CD=bside[ii*atomperresidue+1];
          	BCD=angleside[ii*atomperresidue+1];
          	ABCD=dihside[ii*atomperresidue+1]; 
			//compute HB1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
	   		} 
          	CD=bside[ii*atomperresidue+2];
          	BCD=angleside[ii*atomperresidue+2];
          	ABCD=dihside[ii*atomperresidue+2]; 
			//compute HB2
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
	   		} 
          	CD=bside[ii*atomperresidue+3];
          	BCD=angleside[ii*atomperresidue+3];
          	ABCD=dihside[ii*atomperresidue+3]; 
			//compute CG
          	findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+9]=D[kk];   
	   		} 
			
			//rearrange for computing CD1
         	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}

	   		CD=bside[ii*atomperresidue+4];
          	BCD=angleside[ii*atomperresidue+4];
          	ABCD=dihside[ii*atomperresidue+4]; 
			//compute CD1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+12]=D[kk];   
	   		} 
			
			//rearrange for computing NE1, HD1			
         	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
          	CD=bside[ii*atomperresidue+5];
          	BCD=angleside[ii*atomperresidue+5];
          	ABCD=dihside[ii*atomperresidue+5]; 
			//compute HD1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+15]=D[kk];   
	   		} 
	   		CD=bside[ii*atomperresidue+6];
          	BCD=angleside[ii*atomperresidue+6];
          	ABCD=dihside[ii*atomperresidue+6]; 
			//compute NE1
          	findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+18]=D[kk];   
	   		} 
			
			//rearrange for computing CE2,HE1
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
	   		CD=bside[ii*atomperresidue+7];
          	BCD=angleside[ii*atomperresidue+7];
          	ABCD=dihside[ii*atomperresidue+7]; 
			//compute HE1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+21]=D[kk];   
	   		} 
          	CD=bside[ii*atomperresidue+8];
          	BCD=angleside[ii*atomperresidue+8];
          	ABCD=dihside[ii*atomperresidue+8];
			//compute CE2			
	   		findXYZ_Dihedral(A, B, C, CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+24]=D[kk];   
	   		}  
			
			//rearrange CD2
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
	   		CD=bside[ii*atomperresidue+9];
          	BCD=angleside[ii*atomperresidue+9];
          	ABCD=dihside[ii*atomperresidue+9]; 
			//compute CD2
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
	   		}
			
			//rearrange CE3
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];    
	  		}
	   		CD=bside[ii*atomperresidue+10];
          	BCD=angleside[ii*atomperresidue+10];
          	ABCD=dihside[ii*atomperresidue+10];
			//compute CE3
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
	   		}
			
			//rearrange HE3,CZ3
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk]; 
	 		 }
	   		CD=bside[ii*atomperresidue+11];
          	BCD=angleside[ii*atomperresidue+11];
          	ABCD=dihside[ii*atomperresidue+11]; 
			//compute HE3
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
	   		}
	   		CD=bside[ii*atomperresidue+12];
          	BCD=angleside[ii*atomperresidue+12];
          	ABCD=dihside[ii*atomperresidue+12]; 
			//compute CZ3
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
	   		}
			
			//rearrange HZ3, CH2
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];    
	  		}
	   		CD=bside[ii*atomperresidue+13];
          	BCD=angleside[ii*atomperresidue+13];
          	ABCD=dihside[ii*atomperresidue+13]; 
			//compute HZ3
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
	   		}
	   		CD=bside[ii*atomperresidue+14];
          	BCD=angleside[ii*atomperresidue+14];
          	ABCD=dihside[ii*atomperresidue+14];  
			//compute CH2			
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+42]=D[kk];   
	   		}
			
			//rearrange CZ2, HH2
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];    
	  		}
	   		CD=bside[ii*atomperresidue+15];
          	BCD=angleside[ii*atomperresidue+15];
          	ABCD=dihside[ii*atomperresidue+15]; 
			//compute HH2
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+45]=D[kk];   
	   		}
	   		CD=bside[ii*atomperresidue+16];
          	BCD=angleside[ii*atomperresidue+16];
          	ABCD=dihside[ii*atomperresidue+16]; 
			//compute CZ2
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+48]=D[kk];   
	   		}
						
			//rearrange HZ2
        	for (int kk=0;kk<3;kk++){
          		temp[kk]=A[kk];
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
	   		CD=bside[ii*atomperresidue+17];
          	BCD=angleside[ii*atomperresidue+17];
          	ABCD=dihside[ii*atomperresidue+17]; 
			//compute HZ2
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+51]=D[kk];   
	   		}
		ii++;
		}

		else if (str.compare(PHE)==0){
			cout << PHE<< "\n";
          	for (int kk=0;kk<3;kk++){
          		A[kk]=globalXYZ[(ii+1)*12+3+kk];
          		B[kk]=globalXYZ[(ii+1)*12+kk];
          		D[kk]=globalXYZ[(ii+1)*12+6+kk];
	  		}
	   		AC=bside[ii*atomperresidue];
          	BAC=angleside[ii*atomperresidue];
          	DAC=angleimprside[ii];
          	ABCD=dihside[ii*atomperresidue];
			//CB computation
	   		findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk]=C[kk];   
	  		} 
			
			//rearrange HB1,HB2 computation 
         	for (int kk=0;kk<3;kk++){
          		B[kk]=globalXYZ[(ii+1)*12+3+kk];
          		A[kk]=globalXYZ[(ii+1)*12+kk];   
	  		}
	   		CD=bside[ii*atomperresidue+1];
          	BCD=angleside[ii*atomperresidue+1];
          	ABCD=dihside[ii*atomperresidue+1]; 
			//HB1 computation
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
	   		} 
          	CD=bside[ii*atomperresidue+2];
          	BCD=angleside[ii*atomperresidue+2];
          	ABCD=dihside[ii*atomperresidue+2]; 
			//HB2 computation
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
	   		} 
          	CD=bside[ii*atomperresidue+3];
          	BCD=angleside[ii*atomperresidue+3];
          	ABCD=dihside[ii*atomperresidue+3]; 
			//CG computation
          	findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+9]=D[kk];
	   		} 
			//rearrange for CD1
         	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}

	   		CD=bside[ii*atomperresidue+4];
          	BCD=angleside[ii*atomperresidue+4];
          	ABCD=dihside[ii*atomperresidue+4]; 
			//computation CD1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+12]=D[kk];   
	   		} 
			
			//rearrange for CE1, HD1
         	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
          	CD=bside[ii*atomperresidue+5];
          	BCD=angleside[ii*atomperresidue+5];
          	ABCD=dihside[ii*atomperresidue+5]; 
			//computation HD1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+15]=D[kk];   
	   		} 
	   		CD=bside[ii*atomperresidue+6];
          	BCD=angleside[ii*atomperresidue+6];
          	ABCD=dihside[ii*atomperresidue+6]; 
			//computation CE1
          	findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	   		for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+18]=D[kk];   
	   		} 
			
			//rearrange for CZ, HE1
        	for (int kk=0;kk<3;kk++){
          		A[kk]=B[kk];
          		B[kk]=C[kk];
          		C[kk]=D[kk];   
	  		}
	   		CD=bside[ii*atomperresidue+7];
          	BCD=angleside[ii*atomperresidue+7];
          	ABCD=dihside[ii*atomperresidue+7];
			//computation HE1
	   		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          	for (int kk=0;kk<3;kk++){
          		sideXYZ[ii*atomperresidue*3+kk+21]=D[kk];   
	   		} 
			CD=bside[ii*atomperresidue+8];
			BCD=angleside[ii*atomperresidue+8];
			ABCD=dihside[ii*atomperresidue+8]; 
			//computation CZ
			findXYZ_Dihedral(A, B, C, CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+24]=D[kk];   
			}  
			
			
			//rearrange for HZ, CE2
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+9];
			BCD=angleside[ii*atomperresidue+9];
			ABCD=dihside[ii*atomperresidue+9];
			//computation HZ
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
			}
			CD=bside[ii*atomperresidue+10];
			BCD=angleside[ii*atomperresidue+10];
			ABCD=dihside[ii*atomperresidue+10];
			//computation CE2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
			}
			
			//rearrange for CD2, HE2
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+11];
			BCD=angleside[ii*atomperresidue+11];
			ABCD=dihside[ii*atomperresidue+11];
			//computation HE2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
			}
			CD=bside[ii*atomperresidue+12];
			BCD=angleside[ii*atomperresidue+12];
			ABCD=dihside[ii*atomperresidue+12];
			//computation CD2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
			}
			
			//rearrange for HD2
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+13];
			BCD=angleside[ii*atomperresidue+13];
			ABCD=dihside[ii*atomperresidue+13]; 
			//computation HD2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
			}
			ii++;
		}

		else if (str.compare(TYR)==0){
			cout << TYR<< "\n";
         	for (int kk=0;kk<3;kk++){
				A[kk]=globalXYZ[(ii+1)*12+3+kk];
				B[kk]=globalXYZ[(ii+1)*12+kk];
				D[kk]=globalXYZ[(ii+1)*12+6+kk];
			}
			AC=bside[ii*atomperresidue];
			BAC=angleside[ii*atomperresidue];
			DAC=angleimprside[ii];
			ABCD=dihside[ii*atomperresidue];
			//CB computation 
			findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk]=C[kk];   
			} 
			
			//rearrange for CG, HB1, HB2
			for (int kk=0;kk<3;kk++){
				B[kk]=globalXYZ[(ii+1)*12+3+kk];
				A[kk]=globalXYZ[(ii+1)*12+kk];   
			}
			CD=bside[ii*atomperresidue+1];
			BCD=angleside[ii*atomperresidue+1];
			ABCD=dihside[ii*atomperresidue+1]; 
			//computation HB1
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
			} 
			CD=bside[ii*atomperresidue+2];
			BCD=angleside[ii*atomperresidue+2];
			ABCD=dihside[ii*atomperresidue+2]; 
			//computation HB2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
			} 
			CD=bside[ii*atomperresidue+3];
			BCD=angleside[ii*atomperresidue+3];
			ABCD=dihside[ii*atomperresidue+3]; 
			//computation CG
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+9]=D[kk];   
			} 
			
			//rearrange for CD1
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}

			CD=bside[ii*atomperresidue+4];
			BCD=angleside[ii*atomperresidue+4];
			ABCD=dihside[ii*atomperresidue+4];
			//computation CD1
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+12]=D[kk];   
			} 
			
			//rearrange HD1
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+5];
			BCD=angleside[ii*atomperresidue+5];
			ABCD=dihside[ii*atomperresidue+5];
			//computation HD1
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+15]=D[kk];   
			} 
			CD=bside[ii*atomperresidue+6];
			BCD=angleside[ii*atomperresidue+6];
			ABCD=dihside[ii*atomperresidue+6]; 
			//computation CE1
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+18]=D[kk];   
			} 
			
			//rearrange for HE1, CZ
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+7];
			BCD=angleside[ii*atomperresidue+7];
			ABCD=dihside[ii*atomperresidue+7]; 
			//compute HE1
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+21]=D[kk];   
			} 
			CD=bside[ii*atomperresidue+8];
			BCD=angleside[ii*atomperresidue+8];
			ABCD=dihside[ii*atomperresidue+8]; 
			//compute CZ
			findXYZ_Dihedral(A, B, C, CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+24]=D[kk];   
			}  
			
			//rearrange for OH
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+9];
			BCD=angleside[ii*atomperresidue+9];
			ABCD=dihside[ii*atomperresidue+9];
			//compute OH
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
			}
			//rearrange for compute HH
			for (int kk=0;kk<3;kk++){
				temp[kk]=A[kk];
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+10];
			BCD=angleside[ii*atomperresidue+10];
			ABCD=dihside[ii*atomperresidue+10]; 
			//compute HH
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
			}
			
			//rearrange for compute CE2
			for (int kk=0;kk<3;kk++){
				A[kk]=temp[kk];
				B[kk]=A[kk];
				C[kk]=B[kk];
			}
			CD=bside[ii*atomperresidue+11];
			BCD=angleside[ii*atomperresidue+11];
			ABCD=dihside[ii*atomperresidue+11]; 
			//compute CE2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
			}	
			
			//rearrange for computing HE2
			for (int kk=0;kk<3;kk++){
				temp[kk]=A[kk];
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+12];
			BCD=angleside[ii*atomperresidue+12];
			ABCD=dihside[ii*atomperresidue+12]; 
			//compute HE2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
			}
			CD=bside[ii*atomperresidue+13];
			BCD=angleside[ii*atomperresidue+13];
			ABCD=dihside[ii*atomperresidue+13]; 
			//compute CD2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
			}
			//rearrange HD2 computation 
			for (int kk=0;kk<3;kk++){
				A[kk]=B[kk];
				B[kk]=C[kk];
				C[kk]=D[kk];   
			}
			CD=bside[ii*atomperresidue+14];
			BCD=angleside[ii*atomperresidue+14];
			ABCD=dihside[ii*atomperresidue+14];
			//computing HD2
			findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
			for (int kk=0;kk<3;kk++){
				sideXYZ[ii*atomperresidue*3+kk+42]=D[kk];   
			}
			ii++;
		}


		else  {
			//cout << ii<< "\n";
			for(int i=0;i<26;i++){
       			bond[i]=bside[ii*atomperresidue+i];
       			angle[i]=angleside[ii*atomperresidue+i];
       			dihd[i]=dihside[ii*atomperresidue+i];
     		}

     		for(int k=0;k<3;k++){
       			A[k]=globalXYZ[(ii+1)*12+3+k];
        		B[k]=globalXYZ[(ii+1)*12+k];
        		D[k]=globalXYZ[(ii+1)*12+6+k];	
			    
      		}
     		double DAC=angleimprside[ii];
			//CB calculation using improp dihedral angle
     		findXYZ_ImproperDihedral(A, B, D,bond[0],DAC,angle[0], dihd[0],C); 

			// reorder variables for next round of calculation
      		for(int k=0;k<3;k++){
        		sideXYZ[ii*atomperresidue*3+k]=C[k];
        		temp[k]=A[k];
       			A[k]=B[k];
        		B[k]=temp[k];
      		}
      		// calculating atom 6th of the side chain (CG)
      		findXYZ_Dihedral(A,B,C, bond[6], angle[6], dihd[6],D);
			
     		for(int k=0;k<3;k++){
         		sideXYZ[ii*atomperresidue*3+6*3+k]=D[k];
				temp[k]=D[k];
       		}

			//calculating atom (HB1, HB2, or any atom conected to the CB before going to CG LIKE CG1)
      		for(int i=1;i<3;i++){
          		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          		for(int k=0;k<3;k++){
             		sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
         		}
     		}
			// re-arreng variables for computing secondary connection CB (HG11,HG12,HG13)
      		for(int k=0;k<3;k++){
         		A[k]=B[k];
         		B[k]=C[k];
          		C[k]=D[k];          
      		}
			// calculating secondary connection CB (HG11,HG12,HG13)
			for(int i=3;i<6;i++){
           		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
           		for(int k=0;k<3;k++){
              		sideXYZ[ii*atomperresidue*3+3*i+k]=D[k];
          		}
      		}
			// get the CG as the atom C to go for calculation of next level
       		for(int k=0;k<3;k++){
         		C[k]=temp[k];
       		}
			// calculating atom 12th of the side chain (CD)
      		findXYZ_Dihedral(A,B,C, bond[12], angle[12], dihd[12],D);
      		for(int k=0;k<3;k++){
              	sideXYZ[ii*atomperresidue*3+12*3+k]=D[k];
              	temp[k]=D[k];
      		}
			//calculating atom (HG1, HG2, or any atom conected to the CG before going to CD LIKE CD1)
      		for(int i=8;i>6;i--){
          		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          		for(int k=0;k<3;k++){
              		sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
       			 }
     		}
			// re-arreng variables for calculating secondary connection CG (HG11,HG12,HG13)
      		for(int k=0;k<3;k++){
          		A[k]=B[k];
          		B[k]=C[k];
          		C[k]=D[k];          
      		}
			// calculating secondary connection CG (HG11,HG12,HG13)
      		for(int i=9;i<12;i++){
        		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
         		for(int k=0;k<3;k++){
           			sideXYZ[ii*atomperresidue*3+3*i+k]=D[k];
        			}
     		    }
			// get the CD as the atom C to go for calculation of next level
            for(int k=0;k<3;k++){
         		C[k]=temp[k];
       		}
			//calculating atom (HD1, HD2, CE ) connected to CD
     		for(int i=13;i<16;i++){
          		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          		for(int k=0;k<3;k++){
           			sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
         		}
      		}   
			// re-arreng variables for computing CE connection (HE1,HE2,CZ) going to the next level of side chain structure	
      		for(int k=0;k<3;k++){
          		A[k]=B[k];
          		B[k]=C[k];
         		C[k]=D[k];          
      		}
			//calculating atom (HE1, HE2, CZ) connecting to CE
      		for(int i=16;i<19;i++){
      	  		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
      	  		for(int k=0;k<3;k++){
         			sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
        		}
      		}
			// re-arreng variables for computing atoms connected to CZ 
			for(int k=0;k<3;k++){
          		A[k]=B[k];
          		B[k]=C[k];
          		C[k]=D[k];          
      		}
			//computing atoms connected to CE (HE1, HE2) 
			findXYZ_Dihedral(A,B,C, bond[19], angle[19], dihd[19],D);
      		for(int i=20 ;i>18;i--){
      	  		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
      	  		for(int k=0;k<3;k++){
         			sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
        		}
      		}
			// save the C and re-arreng for computing secondary connection
			for(int k=0;k<3;k++){
              	temp[k]=C[k];
              	C[k]=D[k];
      		}
			
			// computing atoms secondly connected to CE (connected to atom 19th) 
      		for(int i=21;i<23;i++){
        		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
        		for(int k=0;k<3;k++){
         			sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
         	    }
   			}
			
			// retrive C from temp again for computing connected atom to CE (CH) 
      		for(int k=0;k<3;k++){
             	C[k]=temp[k];
      	    }
			// computing CH
      		findXYZ_Dihedral(A,B,C, bond[23], angle[23], dihd[23],D);
      		for(int k=0;k<3;k++){
         		sideXYZ[ii*atomperresidue*3+23*3+k]=D[k];
              	C[k]=D[k];
      		}
			// computing atoms connected to CH 
      		for(int i=24;i<26;i++){
      	  		findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
      	  		for(int k=0;k<3;k++){
         			sideXYZ[ii*atomperresidue*3+i*3+k]=D[k];
         		}
      		}

		}
		ii++;
	}
                                                           
	
	return sideXYZ;
}
        


double *temproryXYZ(bond *batom,angle *angleaatom,dihd *dihatom,int noresidue, double *N_P, double *CA_P, double *C_P){
	double *tempXYZ= new double [noresidue*4*3];
    double *N_temp= new double [3];
    double *C_temp= new double [3];
	double *CA_temp= new double [3];
    double *O_temp= new double [3];
	int coeff = 3*4;	
    for (int idx=0;idx<3;idx++){
        //N xyz
        N_temp[idx] = tempXYZ[idx] = N_P[idx];
	    //CA xyz
	    CA_temp[idx] = tempXYZ[idx+3]= CA_P[idx];
	     //C xyz
	    C_temp[idx] = tempXYZ[idx+6]= C_P[idx];
                 }
    //O xyz
	findXYZ_Dihedral(N_temp,CA_temp,C_temp, batom[0].CO,angleaatom[0].CACO,dihatom[0].NCACO,O_temp);
             	

	for (int i=0;i<noresidue;i++){

        for(int j=0;j<3;j++){
			tempXYZ[(i*coeff)+j]= N_temp[j];
			tempXYZ[(i*coeff)+3+j]= CA_temp[j];
			tempXYZ[(i*coeff)+6+j]= C_temp[j];
			tempXYZ[(i*coeff)+9+j]= O_temp[j];
			    }

		//N xyz
           findXYZ_Dihedral(N_temp,CA_temp,C_temp,batom[i].CN,angleaatom[i].CACN,dihatom[i].NCACN, N_temp);

		//CA xyz
           findXYZ_Dihedral(CA_temp,C_temp,N_temp,batom[i+1].NCA,angleaatom[i].CNCA,dihatom[i].CACNCA, CA_temp);

		//C xyz
           findXYZ_Dihedral(C_temp,N_temp,CA_temp,batom[i+1].CAC,angleaatom[i+1].NCAC,dihatom[i].CNCAC,C_temp );

		//O xyz
		   findXYZ_Dihedral(N_temp,CA_temp,C_temp, batom[i+1].CO,angleaatom[i+1].CACO,dihatom[i+1].NCACO,O_temp);			

	 }
	return tempXYZ;
}
	

	

double  norm(struct XYZ A){
	double nA=sqrt(pow(A.X,2)+pow(A.Y,2)+pow(A.Z,2));
	return nA;

}


XYZ cross(XYZ A,XYZ B){

XYZ CrossResult;
CrossResult.X=A.Y*B.Z-A.Z*B.Y;
CrossResult.Y=-(A.X*B.Z-A.Z*B.X);
CrossResult.Z=A.X*B.Y-A.Y*B.X;

return CrossResult;

}



double *calcrotationmatrix(double BCD,double *N){

	double s=sin(BCD);
	double c=cos(BCD);
	double *rotationmatrix=new double [9];

	rotationmatrix[0]=pow(N[0],2)*(1-c)+c;
	rotationmatrix[1]=N[0]*N[1]*(1-c)-N[2]*s;
	rotationmatrix[2]=N[0]*N[2]*(1-c)+N[1]*s;
	rotationmatrix[3]=N[0]*N[1]*(1-c)+N[2]*s;
	rotationmatrix[4]=pow(N[1],2)*(1-c)+c;
	rotationmatrix[5]=N[1]*N[2]*(1-c)-N[0]*s;
	rotationmatrix[6]=N[0]*N[2]*(1-c)-N[1]*s;
	rotationmatrix[7]=N[2]*N[1]*(1-c)+N[0]*s;
	rotationmatrix[8]=pow(N[2],2)*(1-c)+c;
	
	return rotationmatrix;

}



double *reg_matrixmult(double *A, double *B, int m, int n,int p){
	double *mul_res=new double[m*p];
	for(int i=0;i<m*p;i++)
		mul_res[i]=0;
	
	for(int w=0;w<m;w++){
		for(int c=0;c<p;c++){
			for(int j=0;j<n; j++){
				mul_res[w*p+c]+= A[w*n+j]*B[j*p+c];
			}
		}
	}

		return mul_res;

}


double *reg_cross(double *A,double *B){

	double *CrossResult=new double [3];
	CrossResult[0]=A[1]*B[2]-B[1]*A[2];
	CrossResult[1]=A[2]*B[0]-B[2]*A[0];
	CrossResult[2]=A[0]*B[1]-B[0]*A[1];

	return CrossResult;
}


double *mergrot(double *v1,double *v2,double *v3)
{
   double *rot=new  double [9];
   double x[3],y[3],x1[3],y1[3];

   for(int i=0;i<3;i++){
	   x[i]=v2[i]-v1[i];
	   y[i]=v3[i]-v1[i];
   }
	
    double normx=norm(x);	
   for(int i=0;i<3;i++){
	   x1[i]=x[i]/normx;
   }

   double *x_=Trans(x1,3);
   double *Trans_temp1=reg_matrixmult(x1,x_,3,1,3);
   double *Trans_temp2=reg_matrixmult(Trans_temp1,y,3,3,1);
   
   for(int i=0;i<3;i++){
	   y1[i]=y[i]-Trans_temp2[i];
   }
   double normy=norm(y1);
   for(int i=0;i<3;i++){
	   y1[i]=y1[i]/normy;
	 }

   double *z=reg_cross(x1,y1);
   
   for(int i=0;i<3;i++) {
	   rot[3*i]=x1[i];
	   rot[3*i+1]=y1[i];
	   rot[3*i+2]=z[i];
   }
   return rot;
}


double *matrixmult(double *thetarotation,double *D){

	double *mulD=new double[3];
		mulD[0]=thetarotation[0]*D[0]+thetarotation[1]*D[1]+thetarotation[2]*D[2];
		mulD[1]=thetarotation[3]*D[0]+thetarotation[4]*D[1]+thetarotation[5]*D[2];
		mulD[2]=thetarotation[6]*D[0]+thetarotation[7]*D[1]+thetarotation[8]*D[2];
		return mulD;

}

double  norm(double *A){
	double nA=sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2));
	return nA;

}

double* Trans(double *A, int s){
	double *B= new double [s];
	for(int i=0;i<s;i++){
		B[i]=A[i];
	}

	return B;
}