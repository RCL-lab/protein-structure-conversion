/*####################################################################################################
## FileName:    InternaltoCartesian.cu
## Description: The module implemeted in parallel on GPU and removed the independency, to convert internal to cartesian coordinate.
##				The module reads the intenal coordinate of a protein.
##			    The input files of this program is Bond, Angle, Dihd, imprANGLE_SIDE,BOND_SIDE, ANGLE_SIDE
##				DIHD_SIDE,resname. the output can be two raw data file, backboneXYZ, sideXYZ. 
##				by including ItoC.h header file, we can call writeToPDB.cpp and get the pdb like output. (you only need to uncomment that part)
##				The output is simplepdb.txt.
##			    To run this file just need to compile it by nvcc and then run it:  
## 				nvcc InternaltoCartesian.cu -o InternaltoCartesian_out
##				./InternaltoCartesian_out
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
#include <iomanip> 
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuComplex.h"
#include "ItoC.h"

using namespace std; 
#define pi          3.141592653589793238462643383279502884L

struct bond{
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
};



__host__ __device__ void reg_matrixmult(double *A,double *B,int m,int n,int p,double *C){

 double var;
/*multiplication*/     
	for(int w=0;w<m;w++){
		for(int c=0;c<p;c++){
			var=0;                     
			for(int j=0;j<n; j++){
				var+= A[w*n+j]*B[j*p+c];
			}
			C[w*p+c]=var;
		}
	}
}

__device__ void merge (double *N,double *CA, double *C,double *Translationlocal){
	
	double x[3],x1[3],y[3],y1[3],z[3],temp1[9],temp2[3];
	for(int i=0;i<3;i++){
		x[i]=CA[i]-N[i];
		y[i]=C[i]-N[i];
	}
	double normx=sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
	
	for(int i=0;i<3;i++)
		x1[i]=x[i]/normx;
       
       reg_matrixmult(x1,x1,3,1,3,temp1);
	reg_matrixmult(temp1,y,3,3,1,temp2);

	for(int i=0;i<3;i++)
		y1[i]=y[i]-temp2[i];

       double normy=sqrt(pow(y1[0],2)+pow(y1[1],2)+pow(y1[2],2));

	for(int i=0;i<3;i++)
		y1[i]=y1[i]/normy;
  
	/*cross x1,y1*/
	z[0]=x1[1]*y1[2]-x1[2]*y1[1];
	z[1]=-(x1[0]*y1[2]-x1[2]*y1[0]);
	z[2]=x1[0]*y1[1]-x1[1]*y1[0];
	/*End cross x1,y1*/

	/*Convert 4by4 Translation*/
	   for(int i=0;i<3;i++) {
	   Translationlocal[4*i]=x1[i];
	   Translationlocal[4*i+1]=y1[i];
	   Translationlocal[4*i+2]=z[i];
          	   Translationlocal[4*i+3]=N[i];
	   Translationlocal[12+i]=0;
   	   }	
	   Translationlocal[15]=1;
	
	/*End Convert 4by4 Translation*/

}


__host__ __device__ void findXYZ_Dihedral(double *A, double *B, double *C, double CD, double BCD, double ABCD,double *D){

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

__host__ __device__ double getTorsion (double *A, double *B, double *C, double *D)
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

__host__ __device__ void findXYZ_ImproperDihedralTemp(double *A, double *B, double *D, double AC,double BAC, double ABCD,double *C)
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
__host__ __device__ void findXYZ_ImproperDihedral(double *A, double *B, double *D, double AC, double DAC,double BAC, double ABCD,double *C){
       

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
                if(abs(T-ABCD) < diffABCD)
                {
                     diffABCD=abs(T-ABCD);
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
        
        if(abs(T1-ABCD)<abs(T2-ABCD)){
          for (int i=0;i<3;i++)
              C[i]=C1[i];
       }
        else{
          for (int i=0;i<3;i++)
               C[i]=T1;
          }
        }

}



__global__ void temproryXYZ(struct bond *batom,struct angle *angleatom,struct dihd *dihdatom,double *tempXYZ){
              
              int idx=blockDim.x * blockIdx.x + threadIdx.x;
		double N[3],CA[3],C[3];
              double A_NCAC,B_CAC,B_NCA;

              A_NCAC=angleatom[idx].NCAC;
              B_CAC=batom[idx].CAC;
              B_NCA=batom[idx].NCA;

              //N xyz
              for(int i=0;i<3;i++){
	 N[i]=0;
	tempXYZ[idx*12+i]=N[i];
	}
              
              //CA xyz
              CA[0]=B_NCA;
              CA[1]=0;
              CA[2]=0;
		
              for(int i=0;i<3;i++)
	tempXYZ[idx*12+3+i]=CA[i];

              //C xyz
              C[0]=B_NCA+(B_CAC*cosf(pi-(pi*A_NCAC)/180));
              C[1]=B_CAC*sinf(pi-(pi*A_NCAC)/180);
              C[2]=0;

              for(int i=0;i<3;i++)
	tempXYZ[idx*12+6+i]=C[i];


             //o xyz
              double temp_o[3];
              findXYZ_Dihedral(N,CA,C, batom[idx].CO,angleatom[idx].CACO,dihdatom[idx].NCACO,temp_o);
               
              for(int i=0;i<3;i++)
	tempXYZ[idx*12+9+i]=temp_o[i];  

}

__global__ void TranslationXYZ (struct bond *batom,struct angle *angleaatom,struct dihd *dihatom, double *tempXYZ,int noresidue,double *global_XYZ,double *Translation){

       unsigned int tid = threadIdx.x;
       //unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
       int br,u;
       double CD,BCD,ABCD;
       double CA1[4]; double C1[4]; double N1[4]; double O1[4];
       double N2[4]; double CA2[4]; double C2[4]; double O2[4];

       for (unsigned int m=2; m <=blockDim.x ; m *= 2) {
           if (tid % m == m/2 && tid< noresidue-1 && tid!=0) {
               br=tid+1;
				
               /*New N*/
               CD=batom[br-1].CN;
               BCD=angleaatom[br-1].CACN;
               ABCD=dihatom[br-1].NCACN;
				
               for(int i=0;i<3;i++){
                   N1[i]=tempXYZ[(br-1)*12+i];
                   CA1[i]=tempXYZ[(br-1)*12+3+i];
                   C1[i]=tempXYZ[(br-1)*12+6+i];
	  O1[i]=tempXYZ[br*12+9+i];
	}					
				
               findXYZ_Dihedral(N1,CA1,C1,CD,BCD,ABCD,N2);
					
               for(int i=0;i<3;i++){
	tempXYZ[br*12+i]=N2[i];
				}
			       
              /*End New N*/

              /*New CA*/
              CD=batom[br].NCA;
              BCD=angleaatom[br-1].CNCA;
              ABCD=dihatom[br-1].CACNCA;
              findXYZ_Dihedral(CA1,C1,N2,CD,BCD,ABCD,CA2);
              for(int i=0;i<3;i++)
	tempXYZ[br*12+3+i]=CA2[i];
            /*End New CA*/

            /*New C*/
              CD=batom[br].CAC;
              BCD=angleaatom[br].NCAC;
              ABCD=dihatom[br-1].CNCAC;
              findXYZ_Dihedral(C1,N2,CA2,CD,BCD,ABCD,C2);
              for(int i=0;i<3;i++)
	tempXYZ[br*12+6+i]=C2[i];
            /*End New C*/

            /*New O*/
               double Translationlocal[16];
               merge (N2,CA2,C2,Translationlocal);
               reg_matrixmult(Translationlocal,O1,4,4,1,O2);
               for(int i=0;i<3;i++)
	   tempXYZ[br*12+9+i]=O2[i];
           /*End New O*/
				
           /*Update Translation br*/
             for(int i=0;i<16;i++)
	Translation[br*16+i]=Translationlocal[i];			      
          /*End Update Translation br*/

          /*Update rest of chain Translation*/

              if(m>2){

	u=(br-1)+m/2;
                 for(int i=0;i<3;i++){
                     N1[i]=tempXYZ[u*12+i];
                     CA1[i]=tempXYZ[u*12+3+i];
                     C1[i]=tempXYZ[u*12+6+i];
	    O1[i]=tempXYZ[u*12+9+i];
                  }

                   N1[3]=1;
                   CA1[3]=1;
                   C1[3]=1;
	  O1[3]=1;
					
	   reg_matrixmult(Translationlocal,N1,4,4,1,N2);
	   reg_matrixmult(Translationlocal,CA1,4,4,1,CA2);
	   reg_matrixmult(Translationlocal,C1,4,4,1,C2);
	   reg_matrixmult(Translationlocal,O1,4,4,1,O2);
                       
	   for(int i=0;i<3;i++){
                            tempXYZ[u*12+i]=N2[i];
                            tempXYZ[u*12+3+i]=CA2[i];
                            tempXYZ[u*12+6+i]=C2[i];
	           tempXYZ[u*12+9+i]=O2[i];
	    }
	}
      /*End Find Last in Residue*/
      __syncthreads();	
	}	
     }
	
	
}


__global__ void FinalXYZ(double * Translation, double* tempXYZ){
  unsigned int tid = threadIdx.x;
  double temp[16];
  double N1[4]; double N2[4]; double CA1[4]; double CA2[4]; 
  double C1[4]; double C2[4]; double O1[4]; double O2[4];
  int m=1,u1=0;
  int m_th=ceil(logf(tid)/logf(2));
  int index_trans;

    for(int i=0;i<4;i++){
       if(i==3){
                N1[i]=CA1[i]=C1[i]=O1[i]=1;
       }
      else{
                N1[i]=tempXYZ[tid*12+i];
                CA1[i]=tempXYZ[tid*12+3+i];
                C1[i]=tempXYZ[tid*12+6+i];
                O1[i]=tempXYZ[tid*12+9+i];
       }
}

  for(int k=0;k<=m_th;k++,m=m*2){
     index_trans=(((tid-1)/m)*m)+1;
       if(u1!=index_trans){
          for(int j=0;j<16;j++)
             temp[j]=Translation[index_trans*16+j];
	
             reg_matrixmult(temp,N1,4,4,1,N2);
             reg_matrixmult(temp,CA1,4,4,1,CA2);
             reg_matrixmult(temp,C1,4,4,1,C2);
             reg_matrixmult(temp,O1,4,4,1,O2);
   }
	
       for(int i=0;i<4;i++){
         N1[i]=N2[i];
         CA1[i]=CA2[i];
         C1[i]=C2[i];
         O1[i]=O2[i];
       }
         u1=index_trans;
   }

         for(int i=0;i<3;i++){
            tempXYZ[tid*12+i]=N1[i];
            tempXYZ[tid*12+3+i]=CA1[i];
            tempXYZ[tid*12+6+i]=C1[i];
            tempXYZ[tid*12+9+i]=O1[i];
         }
    


}


__global__ void XYZ_adjust(double* tempXYZ)
{
       double Translation_P[16];
       unsigned int tid = threadIdx.x;
       double atom_1[4],atom_2[4];
       double N_P[3]={18.39299964904785,24.60700035095215,0.7289999723434448};
       double CA_P[3]={19.19499969482422,24.738000869750977,1.9490000009536743};
       double C_P[3]={20.62299919128418,25.194000244140625,1.684000015258789};
	merge(N_P,CA_P,C_P,Translation_P);
       for(int i=0;i<3;i++)
       {
	     atom_1[i]=tempXYZ[tid*3+i];
       }
        atom_1[3]=1;
        reg_matrixmult(Translation_P,atom_1,4,4,1,atom_2);
        for(int i=0;i<3;i++)
       {
	   tempXYZ[tid*3+i]=atom_2[i];
       }
}
__global__ void FinalXYZ_SideChain(double *bside,double *angleimpr,double *angleside,double *dihdside, double* tempXYZ, double *side)
{
     int atomperside=26;
     unsigned int tid = threadIdx.x;
     double bond[26], angle[26], dihd[26];
     double A[3], B[3],C[3],D[3],temp[3];
      

     for(int i=0;i<26;i++)
     {
       bond[i]=bside[tid*atomperside+i];
       angle[i]=angleside[tid*atomperside+i];
       dihd[i]=dihdside[tid*atomperside+i];
     }

     for(int k=0;k<3;k++)
     {
        A[k]=tempXYZ[(tid+1)*12+3+k];
        B[k]=tempXYZ[(tid+1)*12+k];
        D[k]=tempXYZ[(tid+1)*12+6+k];	
      }
     double DAC=angleimpr[tid];

     findXYZ_ImproperDihedral(A, B, D,bond[0],DAC,angle[0], dihd[0],C);
      
      for(int k=0;k<3;k++)
      {
        side[tid*atomperside*3+k]=C[k];
        temp[k]=A[k];
        A[k]=B[k];
        B[k]=temp[k];
      }
      int i;
      findXYZ_Dihedral(A,B,C, bond[6], angle[6], dihd[6],D);
      for(int k=0;k<3;k++)
      {
         side[tid*atomperside*3+6*3+k]=D[k];
         temp[k]=D[k];
       }

      for(i=1;i<3;i++)
      {
          findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          for(int k=0;k<3;k++)
         {
             side[tid*atomperside*3+i*3+k]=D[k];
         }
      }
      for(int k=0;k<3;k++)
      {
          A[k]=B[k];
          B[k]=C[k];
          C[k]=D[k];          
      }

      for(int i=3;i<6;i++)
      {
           findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
           for(int k=0;k<3;k++)
           {
               side[tid*atomperside*3+3*i+k]=D[k];
           }
      }
     
       for(int k=0;k<3;k++)
       {
         	C[k]=temp[k];
       }

      findXYZ_Dihedral(A,B,C, bond[12], angle[12], dihd[12],D);
      for(int k=0;k<3;k++)
      {
              side[tid*atomperside*3+12*3+k]=D[k];
              temp[k]=D[k];
      }

      for(i=8;i>6;i--)
      {
          findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          for(int k=0;k<3;k++)
         {
              side[tid*atomperside*3+i*3+k]=D[k];
        }
      }

      for(int k=0;k<3;k++)
      {
          A[k]=B[k];
          B[k]=C[k];
          C[k]=D[k];          
      }

      for(int i=9;i<12;i++)
      {
         findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
         for(int k=0;k<3;k++)
         {
           side[tid*atomperside*3+3*i+k]=D[k];
         }
      }
            
      for(int k=0;k<3;k++)
      {
         C[k]=temp[k];
       }

      for(i=13;i<16;i++)
      {
          findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
          for(int k=0;k<3;k++)
         {
           side[tid*atomperside*3+i*3+k]=D[k];
         }
      } 
      
      for(int k=0;k<3;k++)
      {
          A[k]=B[k];
          B[k]=C[k];
          C[k]=D[k];          
      }

      for(i=16;i<19;i++)
      {
      	  findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
      	  for(int k=0;k<3;k++)
         {
         	side[tid*atomperside*3+i*3+k]=D[k];
        }
      }

      for(int k=0;k<3;k++)
      {
          A[k]=B[k];
          B[k]=C[k];
          C[k]=D[k];          
      }

      findXYZ_Dihedral(A,B,C, bond[19], angle[19], dihd[19],D);
      for(int k=0;k<3;k++)
      {
              side[tid*atomperside*3+19*3+k]=D[k];
              temp[k]=C[k];
              C[k]=D[k];
      }
      
      for(i=20;i<22;i++)
      {
        findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
        for(int k=0;k<3;k++)
         {
         	side[tid*atomperside*3+i*3+k]=D[k];
         }
      }

      for(int k=0;k<3;k++)
      {
              C[k]=temp[k];
      }
      
      findXYZ_Dihedral(A,B,C, bond[22], angle[22], dihd[22],D);
      for(int k=0;k<3;k++)
      {
         	side[tid*atomperside*3+22*3+k]=D[k];
              C[k]=D[k];
      }
      findXYZ_Dihedral(A,B,C, bond[25], angle[25], dihd[25],D);
      for(int k=0;k<3;k++)
      {
         	side[tid*atomperside*3+25*3+k]=D[k];
      }

      for(i=23;i<25;i++)
      {
      	  findXYZ_Dihedral(A,B,C, bond[i], angle[i], dihd[i],D);
      	  for(int k=0;k<3;k++)
         {
         	side[tid*atomperside*3+i*3+k]=D[k];
         }
      }
}



int main()
{

  ifstream mybond,myangle,mydih,mybond_side,myangle_side,mydih_side,myangle_imprside,name;

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
  int amin=115;
  while (!mybond.eof()){
          getline(mybond,str); 
          line++;
   }

   //int noresidue=line/4;
   int noresidue=128;
   cout<< ceil(log2(float(line/4)))<<"\n";
   if (ceil(log2(float(line/4))) != floor(log2(float(line/4))))
       noresidue=int(pow(2,ceil(log2(float(line/4)))+1));
   cout<<"noresidue:"<<noresidue<<"\n";

   int atomperresidue=26;

   bond *batom=(bond*) malloc(noresidue*sizeof(bond));
   angle *angleatom=(angle*)malloc(noresidue*sizeof(angle));
   dihd *dihatom=(dihd*)malloc(noresidue*sizeof(dihd));
   double *globalXYZ=(double*)malloc(noresidue*12*sizeof(double));
   double *Translation=(double*)malloc(noresidue*16*sizeof(double));

   double *bside=(double*) malloc(amin*atomperresidue*sizeof(double));
   double *angleside=(double*)malloc(amin*atomperresidue*sizeof(double));
   double *angleimprside=(double*)malloc(amin*sizeof(double));
   double *dihside=(double*)malloc(amin*atomperresidue*sizeof(double));
   double *sideXYZ=(double*)malloc(amin*atomperresidue*3*sizeof(double));
   double *checkXYZ=(double*)malloc(amin*atomperresidue*sizeof(double));



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

 if (noresidue!=line/4){
        for(int idx=i; idx<noresidue; idx++){
            batom[idx].NCA = batom[idx].CAC = batom[idx].CO = batom[idx].CN = 0.0;
            angleatom[idx].NCAC = angleatom[idx].CACO = angleatom[idx].CACN = angleatom[idx].CNCA = 0.0;
            dihatom[idx].NCACO = dihatom[idx].NCACN = dihatom[idx].CACNCA = dihatom[idx].CNCAC = 0.0;

         }
}
 


for(int j=i;j<noresidue;j++){
		batom[j].NCA=0;
		batom[j].CAC=0; 
		batom[j].CO=0;
		batom[j].CN=0;


		angleatom[j].NCAC=0; 
		angleatom[j].CACO=0;
		angleatom[j].CACN=0;
		angleatom[j].CNCA=0;

		dihatom[j].NCACO=0;
		dihatom[j].NCACN=0;
		dihatom[j].CACNCA=0;
		dihatom[j].CNCAC=0;
}
   
    // reading the sidechain input
    int i1=0;
    mybond.clear();
    mybond.seekg(0, ios::beg);
    while (i1<amin*26){
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


  mybond.close();
  myangle.close();
  mydih.close();
  mybond_side.close();
  myangle_side.close();
  mydih_side.close();

 
  int u=0;

	for(int j=0;j<noresidue;j++){
	    for(int w=0;w<16;w++){
		u=j*16+w;
		if(u==j*16||u==j*16+5||u==j*16+10||u==j*16+15)
		   Translation[u]=1;
		else
		   Translation[u]=0;
	}
	}


 
  //////////////////////////////////////////////////////
cudaEvent_t start_event, stop_event;

/*copy Mem CPU to GPU*/
   bond  *d_batom;
   angle *d_angleatom;
   dihd  *d_dihatom;
   double *d_tempXYZ;
   double *d_globalXYZ;
   double *d_Translation;
   
   double  *d_bside;
   double  *d_angleside;
   double  *d_angleimprside;
   double  *d_dihside;
   double *d_sideXYZ;
   double *checkbside;

   int   s_block=256;
   float elapsed_time;

   cudaMalloc(&d_batom,s_block*sizeof(bond));
   cudaMalloc(&d_angleatom,s_block*sizeof(angle));
   cudaMalloc(&d_dihatom,s_block*sizeof(dihd));
   cudaMalloc(&d_tempXYZ,12*s_block*sizeof(double));
   cudaMalloc(&d_globalXYZ,12*s_block*sizeof(double));
   cudaMalloc(&d_Translation,16*s_block*sizeof(double));
   
//sidechain copy
   cudaMalloc(&d_bside,amin*atomperresidue*sizeof(double));
   cudaMalloc(&d_angleside,amin*atomperresidue*sizeof(double));
   cudaMalloc(&d_angleimprside,amin*sizeof(double));
   cudaMalloc(&d_dihside,amin*atomperresidue*sizeof(double));
   cudaMalloc(&d_sideXYZ,atomperresidue*amin*3*sizeof(double));
   cudaMalloc(&checkbside,amin*atomperresidue*sizeof(double));


   cudaMemcpy(d_batom, batom, s_block*sizeof(bond), cudaMemcpyHostToDevice);
   cudaMemcpy(d_angleatom, angleatom, s_block*sizeof(angle), cudaMemcpyHostToDevice);
   cudaMemcpy(d_dihatom, dihatom, s_block*sizeof(dihd), cudaMemcpyHostToDevice);
   cudaMemcpy(d_Translation, Translation, s_block*16*sizeof(double), cudaMemcpyHostToDevice);

   //sidechain copy
   cudaMemcpy(d_bside, bside, amin*atomperresidue*sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_angleside, angleside, amin*atomperresidue*sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_dihside, dihside, amin*atomperresidue*sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_angleimprside, angleimprside, amin*sizeof(double), cudaMemcpyHostToDevice);

   dim3 dimBlock(s_block, 1);
   dim3 dimBlock1(s_block*4, 1);
   dim3 dimBlock2(amin,1);
   dim3 dimGrid(1,1);
   cudaEventCreate(&start_event);
   cudaEventCreate(&stop_event);
   cudaEventRecord(start_event, 0);
   temproryXYZ<<<dimGrid, dimBlock>>>(d_batom,d_angleatom,d_dihatom,d_tempXYZ);
   TranslationXYZ<<<dimGrid, dimBlock>>>(d_batom,d_angleatom,d_dihatom,d_tempXYZ,noresidue,d_globalXYZ,d_Translation);
   temproryXYZ<<<dimGrid, dimBlock>>>(d_batom,d_angleatom,d_dihatom,d_tempXYZ);
   FinalXYZ<<<dimGrid,dimBlock>>>(d_Translation, d_tempXYZ);


   XYZ_adjust<<<dimGrid,dimBlock1>>>(d_tempXYZ);
   
//sidechain kernel
   FinalXYZ_SideChain<<<dimGrid,dimBlock2>>>(d_bside,d_angleimprside,d_angleside,d_dihside,d_tempXYZ,d_sideXYZ);

   

cudaEventRecord(stop_event, 0);
cudaEventSynchronize(stop_event);
cudaEventElapsedTime(&elapsed_time, start_event, stop_event);
printf ("Time for the kernel: %f ms\n", elapsed_time);
cudaMemcpy(globalXYZ,d_tempXYZ, s_block*12*sizeof(double), cudaMemcpyDeviceToHost);
cudaMemcpy(sideXYZ,d_sideXYZ, amin*atomperresidue*3*sizeof(double), cudaMemcpyDeviceToHost);


/////	PROLINE & Other Cycle////
   double *A=new double [3];
   double *B=new double [3];
   double *D=new double [3];
   double *C=new double [3];
   double *temp=new double [3];
   double AC,DAC,BAC,ABCD,CD,BCD;
   
   name.open("resname.txt");
   name.seekg(0, ios::beg);
   name.clear();
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
   
   int ii=0;

  if(name.is_open()){
   while (!name.eof()){
       getline(name,str);
       if(str.compare(PRO)==0)
       {
		for (int kk=0;kk<3;kk++){
          A[kk]=globalXYZ[(ii+1)*12+3+kk];
          B[kk]=globalXYZ[(ii+1)*12+kk];
          D[kk]=globalXYZ[(ii+1)*12+6+kk];
	    }
		AC=bside[ii*atomperresidue];
        BAC=angleside[ii*atomperresidue];
        DAC=angleimprside[ii];
        ABCD=dihside[ii*atomperresidue];
		
		// computing CB
		findXYZ_ImproperDihedral(A,B,D, AC, DAC,BAC,ABCD,C);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk]=C[kk];   
		}   
        for (int kk=0;kk<3;kk++){
          B[kk]=globalXYZ[(ii+1)*12+3+kk];
          A[kk]=globalXYZ[(ii+1)*12+kk];   
		}
		CD=bside[ii*atomperresidue+1];
        BCD=angleside[ii*atomperresidue+1];
        ABCD=dihside[ii*atomperresidue+1];
		
		//computing CG
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
		} 
		CD=bside[ii*atomperresidue+2];
        BCD=angleside[ii*atomperresidue+2];
        ABCD=dihside[ii*atomperresidue+2];
		//computing CD
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
		for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
		} 
		
		//ADD HB,HG,HD
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
    }
       
       if(str.compare(HIS)==0)
       {
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
        }   
        for (int kk=0;kk<3;kk++){
          B[kk]=globalXYZ[(ii+1)*12+3+kk];
          A[kk]=globalXYZ[(ii+1)*12+kk];   
        }
        CD=bside[ii*atomperresidue+1];
        BCD=angleside[ii*atomperresidue+1];
        ABCD=dihside[ii*atomperresidue+1];          
        findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
	    } 
        CD=bside[ii*atomperresidue+2];
        BCD=angleside[ii*atomperresidue+2];
        ABCD=dihside[ii*atomperresidue+2]; 
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
		for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
		} 
        CD=bside[ii*atomperresidue+3];
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
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
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
		} 
		//rearrenge for calculation CD2		
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
		}
		CD=bside[ii*atomperresidue+9];
        BCD=angleside[ii*atomperresidue+9];
        ABCD=dihside[ii*atomperresidue+9];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
		}
		
		
		//rearrenge for calculation HD2
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
       }
 
       if(str.compare(PHE)==0)
       {
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
	    }   
        for (int kk=0;kk<3;kk++){
          B[kk]=globalXYZ[(ii+1)*12+3+kk];
          A[kk]=globalXYZ[(ii+1)*12+kk];   
	    }
	    CD=bside[ii*atomperresidue+1];
        BCD=angleside[ii*atomperresidue+1];
        ABCD=dihside[ii*atomperresidue+1];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
        for (int kk=0;kk<3;kk++){
        sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
	    } 
        CD=bside[ii*atomperresidue+2];
        BCD=angleside[ii*atomperresidue+2];
        ABCD=dihside[ii*atomperresidue+2]; 
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
	    for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
	    } 
        CD=bside[ii*atomperresidue+3];
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
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
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
	    }  
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+9];
        BCD=angleside[ii*atomperresidue+9];
        ABCD=dihside[ii*atomperresidue+9];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+10];
        BCD=angleside[ii*atomperresidue+10];
        ABCD=dihside[ii*atomperresidue+10];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+11];
        BCD=angleside[ii*atomperresidue+11];
        ABCD=dihside[ii*atomperresidue+11];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+12];
        BCD=angleside[ii*atomperresidue+12];
        ABCD=dihside[ii*atomperresidue+12];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+13];
        BCD=angleside[ii*atomperresidue+13];
        ABCD=dihside[ii*atomperresidue+13];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
	    }
      }

       if(str.compare(TYR)==0)
       {
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
		}   
        for (int kk=0;kk<3;kk++){
          B[kk]=globalXYZ[(ii+1)*12+3+kk];
          A[kk]=globalXYZ[(ii+1)*12+kk];   
		}
		CD=bside[ii*atomperresidue+1];
        BCD=angleside[ii*atomperresidue+1];
        ABCD=dihside[ii*atomperresidue+1];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
		} 
        CD=bside[ii*atomperresidue+2];
        BCD=angleside[ii*atomperresidue+2];
        ABCD=dihside[ii*atomperresidue+2]; 
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
		for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
		}	 
        CD=bside[ii*atomperresidue+3];
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
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
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
		}  
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
		}
		CD=bside[ii*atomperresidue+9];
        BCD=angleside[ii*atomperresidue+9];
        ABCD=dihside[ii*atomperresidue+9];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
		}
        for (int kk=0;kk<3;kk++){
          temp[kk]=A[kk];
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
		}
		CD=bside[ii*atomperresidue+10];
        BCD=angleside[ii*atomperresidue+10];
        ABCD=dihside[ii*atomperresidue+10];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
		}
        for (int kk=0;kk<3;kk++){
          A[kk]=temp[kk];
          B[kk]=A[kk];
          C[kk]=B[kk];
		}
		CD=bside[ii*atomperresidue+11];
        BCD=angleside[ii*atomperresidue+11];
        ABCD=dihside[ii*atomperresidue+11];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
        sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
	   }
        for (int kk=0;kk<3;kk++){
          temp[kk]=A[kk];
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
		}
		CD=bside[ii*atomperresidue+12];
        BCD=angleside[ii*atomperresidue+12];
        ABCD=dihside[ii*atomperresidue+12];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
		}
		CD=bside[ii*atomperresidue+13];
        BCD=angleside[ii*atomperresidue+13];
        ABCD=dihside[ii*atomperresidue+13];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
		}
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
		}
		CD=bside[ii*atomperresidue+14];
        BCD=angleside[ii*atomperresidue+14];
        ABCD=dihside[ii*atomperresidue+14];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
		for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+42]=D[kk];   
		}
       }

       if(str.compare(TRP)==0)
       {
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
		}   
        for (int kk=0;kk<3;kk++){
          B[kk]=globalXYZ[(ii+1)*12+3+kk];
          A[kk]=globalXYZ[(ii+1)*12+kk];   
		}
		CD=bside[ii*atomperresidue+1];
        BCD=angleside[ii*atomperresidue+1];
        ABCD=dihside[ii*atomperresidue+1];          
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
          
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+3]=D[kk];   
		} 
        CD=bside[ii*atomperresidue+2];
        BCD=angleside[ii*atomperresidue+2];
        ABCD=dihside[ii*atomperresidue+2]; 
		findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
		for (int kk=0;kk<3;kk++){
        sideXYZ[ii*atomperresidue*3+kk+6]=D[kk];   
		} 
        CD=bside[ii*atomperresidue+3];
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
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
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
	    }  
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+9];
        BCD=angleside[ii*atomperresidue+9];
        ABCD=dihside[ii*atomperresidue+9];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+27]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];    
	   }
	    CD=bside[ii*atomperresidue+10];
        BCD=angleside[ii*atomperresidue+10];
        ABCD=dihside[ii*atomperresidue+10];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+30]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk]; 
	   }
	    CD=bside[ii*atomperresidue+11];
        BCD=angleside[ii*atomperresidue+11];
        ABCD=dihside[ii*atomperresidue+11];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+33]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+12];
        BCD=angleside[ii*atomperresidue+12];
        ABCD=dihside[ii*atomperresidue+12];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+36]=D[kk];   
	   }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];    
	   }
	    CD=bside[ii*atomperresidue+13];
        BCD=angleside[ii*atomperresidue+13];
        ABCD=dihside[ii*atomperresidue+13];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+39]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+14];
        BCD=angleside[ii*atomperresidue+14];
        ABCD=dihside[ii*atomperresidue+14];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+42]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];    
	    }
	    CD=bside[ii*atomperresidue+15];
        BCD=angleside[ii*atomperresidue+15];
        ABCD=dihside[ii*atomperresidue+15];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+45]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+16];
        BCD=angleside[ii*atomperresidue+16];
        ABCD=dihside[ii*atomperresidue+16];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+48]=D[kk];   
	    }
        for (int kk=0;kk<3;kk++){
          temp[kk]=A[kk];
          A[kk]=B[kk];
          B[kk]=C[kk];
          C[kk]=D[kk];   
	    }
	    CD=bside[ii*atomperresidue+17];
        BCD=angleside[ii*atomperresidue+17];
        ABCD=dihside[ii*atomperresidue+17];          
	    findXYZ_Dihedral(A,B,C,CD,BCD,ABCD,D);
        for (int kk=0;kk<3;kk++){
          sideXYZ[ii*atomperresidue*3+kk+51]=D[kk];   
	    }
        
       }

       
       ii++;
   }
  }
  
  name.close();


 
 
  cudaFree(d_tempXYZ);
   
   
   cudaFree(d_batom);
   cudaFree(d_angleatom);
   cudaFree(d_dihatom);
   cudaFree(d_Translation);
   cudaFree(d_globalXYZ);
  cudaFree(d_bside);
  cudaFree(d_angleside);
  cudaFree(d_sideXYZ);
  cudaFree(d_dihside);


   ofstream myfile ("finalXYZ.txt");
   if (myfile.is_open())
   {
	   for(int k=0;k<noresidue;k++){

		   myfile<<"N, "<<globalXYZ[k*12+0]<<","<<globalXYZ[k*12+1]<<","<<globalXYZ[k*12+2]<<"\n";
		   myfile<<"CA, "<<globalXYZ[k*12+3]<<","<<globalXYZ[k*12+4]<<","<<globalXYZ[k*12+5]<<"\n";
		   myfile<<"C, "<<globalXYZ[k*12+6]<<","<<globalXYZ[k*12+7]<<","<<globalXYZ[k*12+8]<<"\n";
		   myfile<<"O,"<<globalXYZ[k*12+9]<<","<<globalXYZ[k*12+10]<<","<<globalXYZ[k*12+11]<<"\n";

	   }
    myfile.close();
   }
  else cout << "Unable to open file";
 

//sidechain output
   ofstream myfile_s ("side_XYZ.txt");
   if (myfile_s.is_open())
   {
	   for(int k=0;k<amin;k++){
              for(int u=0;u<26;u++){
		    myfile_s<<u<<" "<<sideXYZ[k*atomperresidue*3+3*u]<<","<<sideXYZ[k*atomperresidue*3+3*u+1]<<","<<sideXYZ[k*atomperresidue*3+3*u+2]<<"\n";
	   }
      }
    myfile_s.close();
   }
  else cout << "Unable to open file";

  return 0

}


