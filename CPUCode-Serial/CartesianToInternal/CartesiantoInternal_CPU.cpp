/*####################################################################################################
## FileName:    CartesiantoInternal_CPU.cpp
## Description: The module read the Cartesian coordinate of a protein and the list of atoms that make bond, angle and dihedral together. 
##				serially convert it to Internal coordinate by using formula.
##				The input files of this program is trajectory.txt, bond.txt, angle.txt, improper.txt, proper.txt
##				The output is 4 different output file for each internal Cartesian coordinate.
##			    To run this file just need to compile it by g++ and then run it:  
## 				g++ CartesiantoInternal_CPU.cpp -o CartesiantoInternal_out
##				./CartesiantoInternal_out
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
#include <vector>                       
#include <sys/time.h>
#include <fstream>
#include <cstring>
#include <string>
#include <istream>
#include <math.h>

using namespace std;


#define Natom  2443    
#define Nframe 1000    

double *BondDistance,*newTheta,*Torsion1,*Torsion2;


struct Trajectory {
	double x;
    double y;
	double z;
};

struct timeval start, end;
long mtime, secondtraj, usecondtraj,secondb, usecondb, seconda,useconda, secondd,usecondd, secondim,usecondim,seconds,useconds;
void populateVector ( std:: vector< int > &A, double m );             
void printVector (const std:: vector< int > A ); 
Trajectory *readTrajectory(double *data);
void getBondDistance(Trajectory *traj,int *bonds, int number_bonds, double *BondDistance);
void getAngle(Trajectory *traj,int *angles,int number_angles, double *newTheta);
void getTorsion(Trajectory *traj,int *dihedrals, int number_dihedrals, double *Torsion);

int main()
{





	vector<int> atomlist (Natom);
	vector<int>::iterator i;

		for (i=atomlist.begin(); i!= atomlist.end();i++)
			 {
				 atomlist[i-atomlist.begin()]=i+1-atomlist.begin();
			 }
	double *data;
//get the number of lines in file to get the number of bonds, angles and dihedrals

int number_bonds = 0,number_angles=0,number_dihedrals=0,number_impropers=0,number_traj=0,b=0,a=0,ip=0,d=0,t=0;




ifstream mybond("bond.txt");
ifstream myangle("theta.txt");
ifstream mydihedral("propers.txt");
ifstream myimproper("impropers.txt");
ifstream mytraj("trajectory.txt");


string line;
//***** BONDLENGTH file

	while (!mybond.eof())
	{
		getline (mybond,line);
		number_bonds++;
	}

	cout << "Number of lines in text file: " << number_bonds<<"\n";
	int *bonds=(int*) malloc (number_bonds*sizeof(int));
	mybond.clear();
	mybond.seekg(0, ios::beg);

	while (!mybond.eof())
	{

		getline (mybond,line);
		bonds[b]=stoi(line);
		b++;
	}


//*** Angles file

	while (!myangle.eof() )
	{
		getline (myangle,line);
		number_angles++;
	}
	cout << "Number of lines in text file: " << number_angles<<"\n";
	int *angles=(int*) malloc (number_angles*sizeof(int));
	myangle.clear();
	myangle.seekg(0, ios::beg);

	while (!myangle.eof() )
	{
		getline (myangle,line);
		angles[a]=stoi(line);
		a++;
	}



//dihedral file

	while (!mydihedral.eof() )
	{
		getline (mydihedral,line);
		number_dihedrals++;
	}

	cout << "Number of lines in text file: " << number_dihedrals<<"\n";

	int *dihedrals=(int*) malloc (number_dihedrals*sizeof(int));

	mydihedral.clear();
	mydihedral.seekg(0, ios::beg);

	while (!mydihedral.eof() )
	{
		getline (mydihedral,line);
		dihedrals[d]=stoi(line);
		d++;
	}


//improper file

	while (!myimproper.eof() )
	{
		getline (myimproper,line);
		number_impropers++;
	}

	cout << "Number of lines in text file: " << number_impropers<<"\n";
	int *impropers=(int*) malloc (number_impropers*sizeof(int));

	myimproper.clear();
	myimproper.seekg(0, ios::beg);

	while (!myimproper.eof() )
	{
		getline (myimproper,line);
		impropers[ip]=stoi(line);
		ip++;
	}


	


//trajectory file

	while (!mytraj.eof() )
	{
		getline (mytraj,line);
		number_traj++;
	}

	cout << "Number of lines in text file: " << number_traj<<"\n";
	data=(double*) malloc (number_traj*sizeof(double));
	if(data!=NULL)  printf("correct");
	mytraj.clear();
	mytraj.seekg(0, ios::beg);
	while (!mytraj.eof())
	{
		getline (mytraj,line);
		data[t]=stof(line);
            //  cout<<t<<"\n";
		t++;
	}


//////*******call functions
 
struct Trajectory *traj=(Trajectory*) malloc (Natom*Nframe*3*sizeof(double));
//gettimeofday(&start, NULL);
traj=readTrajectory(data);
//gettimeofday(&end, NULL);
 //secondtraj  = end.tv_sec  - start.tv_sec;
 //usecondtraj = end.tv_usec - start.tv_usec;
 //mtime = ((secondtraj) * 1000 + usecondtraj/1000.0) + 0.5;
//printf("Elapsed time traj: %ld milliseconds\n", mtime);

gettimeofday(&start, NULL);
getBondDistance(traj,bonds,number_bonds,BondDistance);
gettimeofday(&end, NULL);
cout<<number_bonds/2<<"\n";
secondb  = end.tv_sec  - start.tv_sec;
usecondb = end.tv_usec - start.tv_usec;
mtime = ((secondb) * 1000 + usecondb/1000.0) + 0.5;
secondb=start.tv_sec;
usecondb=start.tv_usec;
printf("Elapsed time bond: %ld milliseconds\n", mtime);

gettimeofday(&start, NULL);
getAngle(traj,angles,number_angles,newTheta);
gettimeofday(&end, NULL);
cout<<number_angles/3<<"\n";
seconda  = end.tv_sec  - start.tv_sec;
useconda = end.tv_usec - start.tv_usec;
mtime = ((seconda) * 1000 + useconda/1000.0) + 0.5;
printf("Elapsed time angle: %ld milliseconds\n", mtime);

gettimeofday(&start, NULL);
getTorsion(traj,dihedrals,number_dihedrals, Torsion1);
gettimeofday(&end, NULL);
cout<<number_dihedrals/4<<"\n";
secondd  = end.tv_sec  - start.tv_sec;
usecondd = end.tv_usec - start.tv_usec;
mtime = ((secondd) * 1000 + usecondd/1000.0) + 0.5;
printf("Elapsed time dihed: %ld milliseconds\n", mtime);

gettimeofday(&start, NULL);
getTorsion(traj,impropers,number_impropers, Torsion2);
gettimeofday(&end, NULL);
cout<<number_impropers/4<<"\n";
 secondim  = end.tv_sec  - start.tv_sec;
 usecondim = end.tv_usec - start.tv_usec;
 mtime = ((secondim) * 1000 + usecondim/1000.0) + 0.5;
printf("Elapsed time improp: %ld milliseconds\n", mtime);
secondim=end.tv_sec;
usecondim = end.tv_usec;


//
///**write the result back/////


//////////time calculation////

  /*  seconds  = secondim-secondb;
    useconds = usecondim -usecondb;

    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    printf("Elapsed time: %ld milliseconds\n", mtime);*/
///////***********************
	//system("pause");
		return 0;


}

Trajectory *readTrajectory(double *data)
{
	
	struct Trajectory *traj=(Trajectory*) malloc (Natom*Nframe*3*sizeof(double));

		for(int i=0;i<Natom;i++){
			for(int j=0;j<Nframe;j++){

			traj[j*Natom+i].x=data[j*3*Natom+(i*3)];
			traj[j*Natom+i].y=data[j*3*Natom+(i*3+1)];
			traj[j*Natom+i].z=data[j*3*Natom+(i*3+2)];

	     	}
		  }
		return traj;


}

void getBondDistance(Trajectory *traj,int *bonds, int number_bonds ,double *BondDistance )
{
	BondDistance=(double*) malloc ((number_bonds/2)*Nframe*sizeof(double));
	int A,B;

	for(int i=0;i<Nframe;i++){
		for(int j=0;j<number_bonds;j=j+2){

		A=bonds[j];
		B=bonds[j+1];
		double disx=pow(traj[i*Natom+(A-1)].x-traj[i*Natom+(B-1)].x,2);
		double disy=pow(traj[i*Natom+(A-1)].y-traj[i*Natom+(B-1)].y,2);
		double disz=pow(traj[i*Natom+(A-1)].z-traj[i*Natom+(B-1)].z,2);
		BondDistance[i*(number_bonds/2)+(j/2)]=sqrt(disx+disy+disz);
		}
	
	}

/*FILE *Fb=fopen("/home/scratch/mahsa/bondd.txt","w");
for(int h=0;h<(number_bonds/2)*Nframe;h++) fprintf(Fb, "%f\n",BondDistance[h]);
fclose(Fb);*/
}



void getAngle(Trajectory *traj,int *angles, int number_angles, double *newTheta){
		

		double VBA[3],VBC[3],dotprod,magBA,magBC;//x,y,z multiply by 3 
		int A,B,C;
		
		

		//searching for sets of angles to be in the range of atom no.

		newTheta=(double*) malloc ((number_angles/3)*Nframe*sizeof(double));
		////calculating Angles 
		
		for(int j=0;j<number_angles;j=j+3){
		 for(int i=0;i<Nframe;i++){

			A=angles[j];
			B=angles[j+1];
			C=angles[j+2];

			VBA[0]=traj[i*Natom+(B-1)].x-traj[i*Natom+(A-1)].x;
			VBA[1]=traj[i*Natom+(B-1)].y-traj[i*Natom+(A-1)].y;
			VBA[2]=traj[i*Natom+(B-1)].z-traj[i*Natom+(A-1)].z;

			

			VBC[0]=traj[i*Natom+(B-1)].x-traj[i*Natom+(C-1)].x;
			VBC[1]=traj[i*Natom+(B-1)].y-traj[i*Natom+(C-1)].y;
			VBC[2]=traj[i*Natom+(B-1)].z-traj[i*Natom+(C-1)].z;


			dotprod=VBA[0]*VBC[0]+VBA[1]*VBC[1]+VBA[2]*VBC[2];
			magBA=sqrt(pow(VBA[0],2)+pow(VBA[1],2)+pow(VBA[2],2));
			magBC=sqrt(pow(VBC[0],2)+pow(VBC[1],2)+pow(VBC[2],2));
			newTheta[i*(number_angles/3)+j/3]=acos(dotprod/(magBA*magBC));

	
		}
	}

//gettimeofday(&end, NULL);

/*FILE *Fa=fopen("home/scratch/mahsa/ang.txt","w");
for(int h=0;h<(number_angles/3)*Nframe;h++) fprintf(Fa, "%f\n",newTheta[h]);
fclose(Fa);*/
}



void getTorsion(Trajectory *traj,int *dihedrals,int number_dihedrals, double *Torsion)
{
		
		int A,B,C,D;
		double VAB[3],VBC[3],VCD[3],norm,term1[3],term2[3],term3[3],dotpro1,dotpro2;
		
		Torsion=(double*) malloc ((number_dihedrals/4)*Nframe*sizeof(double));
		//calculating the Torsion

		for(int j=0;j<number_dihedrals;j=j+4){

			A=dihedrals[j];
		    B=dihedrals[j+1];
			C=dihedrals[j+2];
			D=dihedrals[j+3];

		for(int i=0;i<Nframe;i++){
		

			VAB[0]=traj[i*Natom+(B-1)].x-traj[i*Natom+(A-1)].x;
			VAB[1]=traj[i*Natom+(B-1)].y-traj[i*Natom+(A-1)].y;
			VAB[2]=traj[i*Natom+(B-1)].z-traj[i*Natom+(A-1)].z;

			VBC[0]=traj[i*Natom+(C-1)].x-traj[i*Natom+(B-1)].x;
			VBC[1]=traj[i*Natom+(C-1)].y-traj[i*Natom+(B-1)].y;
			VBC[2]=traj[i*Natom+(C-1)].z-traj[i*Natom+(B-1)].z;

			VCD[0]=traj[i*Natom+(D-1)].x-traj[i*Natom+(C-1)].x;
			VCD[1]=traj[i*Natom+(D-1)].y-traj[i*Natom+(C-1)].y;
			VCD[2]=traj[i*Natom+(D-1)].z-traj[i*Natom+(C-1)].z;
             
			norm=sqrt(pow(VBC[0],2)+pow(VBC[1],2)+pow(VBC[2],2));
			term1[0]=norm*VAB[0];
			term1[1]=norm*VAB[1];
			term1[2]=norm*VAB[2];

			term2[0]=(VBC[1]*VCD[2])-(VBC[2]*VCD[1]);
			term2[1]=-((VBC[0]*VCD[2])-(VBC[2]*VCD[0]));
			term2[2]=(VBC[0]*VCD[1])-(VBC[1]*VCD[0]);

			dotpro1=term1[0]*term2[0]+term1[1]*term2[1]+term1[2]*term2[2];

			term3[0]=(VAB[1]*VBC[2])-(VAB[2]*VBC[1]);
			term3[1]=-((VAB[0]*VBC[2])-(VAB[2]*VBC[0]));
			term3[2]=(VAB[0]*VBC[1])-(VAB[1]*VBC[0]);
		
			dotpro2=term2[0]*term3[0]+term2[1]*term3[1]+term2[2]*term3[2];
			Torsion[i*(number_dihedrals/4)+(j/4)]=atan2(dotpro1,dotpro2);
			
		}
	}




}


void populateVector (vector< int > &A, double m )

	{	
		 
		 int ranum;
		 time_t second; 
		 srand( time(&second));

		 vector<int>::iterator i;

		 for (i=A.begin(); i!= A.end();i++)
			 {
				 ranum= rand()% int(m);
				 A[i-A.begin()]=ranum;
			 }
 
	}


void printVector(vector< int >A)
	{

		vector<int>::iterator ii;
		
		for( ii = A.begin(); ii !=A.end(); ii++)
			{
				cout << *ii <<";"<<"\n";
			}
				cout << endl<<endl;

	}