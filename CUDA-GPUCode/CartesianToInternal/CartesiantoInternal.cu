#include "headerReadTrajectory.h"
#include <cuda.h>
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

static void HandleError( cudaError_t err, const char *file, int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}


__global__ void Bdistance(int* d_bond,float* d_data,int pd_bonds,float* d_BondDistance){
	int idx=blockDim.x * blockIdx.x + threadIdx.x;
	int idy=blockDim.y * blockIdx.y + threadIdx.y;
	float dis=0;

	int A=d_bond[idx*2];
	int B=d_bond[idx*2+1];

       int iAdata=idy*Natom*3+(A-1)*3;
	int iBdata=idy*Natom*3+(B-1)*3;

	float disx=d_data[iAdata]-d_data[iBdata];
	float disy=d_data[iAdata+1]-d_data[iBdata+1];
	float disz=d_data[iAdata+2]-d_data[iBdata+2];
	
	disx=pow(disx,2);
	disy=pow(disy,2);
	disz=pow(disz,2);
	
	dis=sqrt(disy+disx+disz);
	
	d_BondDistance[idy*pd_bonds+idx]=dis;

	
}

__global__ void Angle(int* d_angles,float* d_data,int pd_angles,float* d_newTheta){
	int idx=blockDim.x * blockIdx.x + threadIdx.x;
	int idy=blockDim.y * blockIdx.y + threadIdx.y;

	
	int A=d_angles[idx*3];
	int B=d_angles[idx*3+1];
	int C=d_angles[idx*3+2];
	float VBA[3],VBC[3],dotprod=0,magBA=0,magBC=0,data_tempB[3];

		
      	

	//#pragma unroll
	for(int i=0;i<3;i++){
///Access mem once
	data_tempB[i]=d_data[idy*Natom*3+(B-1)*3+i];

	VBA[i]=data_tempB[i]-d_data[idy*Natom*3+(A-1)*3+i];
	VBC[i]=data_tempB[i]-d_data[idy*Natom*3+(C-1)*3+i];
       dotprod+=VBA[i]*VBC[i];
	magBA+=pow(VBA[i],2);
	magBC+=pow(VBC[i],2);
	 }

	magBA=sqrt(magBA);
	magBC=sqrt(magBC);
	d_newTheta[idy*pd_angles+idx]=acos(dotprod/(magBA*magBC));
	
	

}

__global__ void Torsion(int* d_dihedrals,float* d_data,int pd_dihedrals,float* d_Torsion){
	int idx=blockDim.x * blockIdx.x + threadIdx.x;
	int idy=blockDim.y * blockIdx.y + threadIdx.y;
	

	int A=d_dihedrals[idx*4];
	int B=d_dihedrals[idx*4+1];
	int C=d_dihedrals[idx*4+2];
	int D=d_dihedrals[idx*4+3];
	float VAB[3],VBC[3],VCD[3],term1[3],term2[3],term3[3],norm=0,dotpro1,dotpro2;
	float tempBC[6];
	

       tempBC[0]=d_data[idy*Natom*3+(B-1)*3];
	tempBC[1]=d_data[idy*Natom*3+(B-1)*3+1];
	tempBC[2]=d_data[idy*Natom*3+(B-1)*3+2];

       tempBC[3]=d_data[idy*Natom*3+(C-1)*3];
	tempBC[4]=d_data[idy*Natom*3+(C-1)*3+1];
	tempBC[5]=d_data[idy*Natom*3+(C-1)*3+2];

       VAB[0]=tempBC[0]-d_data[idy*Natom*3+(A-1)*3];
	VAB[1]=tempBC[1]-d_data[idy*Natom*3+(A-1)*3+1];
	VAB[2]=tempBC[2]-d_data[idy*Natom*3+(A-1)*3+2];

	VBC[0]=tempBC[3]-tempBC[0];
	VBC[1]=tempBC[4]-tempBC[1];
	VBC[2]=tempBC[5]-tempBC[2];
	
	VCD[0]=d_data[idy*Natom*3+(D-1)*3]-tempBC[3];
	VCD[1]=d_data[idy*Natom*3+(D-1)*3+1]-tempBC[4];
	VCD[2]=d_data[idy*Natom*3+(D-1)*3+2]-tempBC[5];




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

	d_Torsion[idy*pd_dihedrals+idx]=atan2(dotpro1,dotpro2);
	
	

}




int main()
{


	

	vector<int> atomlist (Natom);
	vector<int>::iterator i;

		for (i=atomlist.begin(); i!= atomlist.end();i++)
			 {
				 atomlist[i-atomlist.begin()]=i-atomlist.begin();
			 }

//get the number of lines in file to get the number of bonds, angles and dihedrals

int number_bonds=0,number_angles=0,number_dihedrals=0,number_impropers=0,number_traj=0;
string line;
//***** BONDLENGTH file
number_bonds=getLinebondfile();

int pd_bonds=ceil(log2(float(number_bonds/2)));
pd_bonds=pow(2,pd_bonds);
int *bonds=(int*) malloc (pd_bonds*2*sizeof(int));

rfbond(pd_bonds,bonds);


//*** Angles file
number_angles=getLineanglefile();

int pd_angles=ceil(log2(float(number_angles/3)));
pd_angles=pow(2,pd_angles);
int *angles=(int*) malloc (pd_angles*3*sizeof(int));

rfangle(pd_angles,angles);

//dihedral file
number_dihedrals=getLineproperfile();

int pd_dihedrals=ceil(log2(float(number_dihedrals/4)));
pd_dihedrals=pow(2,pd_dihedrals);
int *dihedrals=(int*) malloc (pd_dihedrals*4*sizeof(int));

rfproper(pd_dihedrals,dihedrals);


//improper file
number_impropers=getLineimproperfile();

int pd_impropers=ceil(log2(float(number_impropers/4)));
pd_impropers=pow(2,pd_impropers);
int *impropers=(int*) malloc (pd_impropers*4*sizeof(int));

rfimproper(pd_impropers,impropers);

//trajectory file
number_traj=getLinetrajfile();

int pd_frame=ceil(log2(float(Nframe)));
pd_frame=pow(2,pd_frame);
float *data=(float*) malloc (Natom*pd_frame*3*sizeof(float));

rftrajectory(pd_frame,data);


//////*******padding
cout<<"pd_bonds:"<<pd_bonds<<"\npd_angles"<<pd_angles;
cout<<"\npd_propers"<<pd_dihedrals<<"\npd_impropers:"<<pd_impropers<<"\npd_frame"<<pd_frame;

/////CUDA PART////////////////

int *d_bond,*d_angle,*d_dih,*d_impro;
float *d_data,*d_BondDistance,*d_newTheta,*d_Torsion1,*d_Torsion2;

int FullData=pd_bonds+pd_angles+pd_dihedrals+pd_impropers+pd_bonds*pd_frame+pd_angles*pd_frame+pd_dihedrals*pd_frame+pd_impropers*pd_frame;
int chunk=0;
int seg_bonds=256, seg_angles=256,seg_dihedrals=256,seg_impropers=256, seg_frame=256;
int offset_bonds=0, offset_angles=0, offset_dihedrals=0, offset_Distance=0, 
offset_Theta=0, offset_Torsion1=0,offset_impropers=0, offset_Torsion2=0;



int TB=32;

dim3 dimBlockB(TB,TB),dimBlockA(TB,TB),dimBlockP(TB,TB),dimBlockIP(4,TB);

dim3 dimGridB(seg_bonds/TB,seg_frame/TB),dimGridA(seg_angles/TB,seg_frame/TB),
     dimGridP(seg_dihedrals/TB,seg_frame/TB),dimGridIP(seg_impropers/TB,seg_frame/TB);
/////////////allocate Memory////////////
cudaError_t err;
int deviceCount;

cudaDeviceReset();
cudaGetDeviceCount(&deviceCount);
cout<<"device:"<<deviceCount<<"\n";
cudaDeviceProp deviceProp ;
cudaGetDeviceProperties ( &deviceProp , 0) ;

cout << "Using CUDA device" << deviceProp.name << endl ;

cudaEvent_t start, stop;
float time;
cudaEventCreate(&start);
cudaEventCreate(&stop);


cudaMalloc(&d_bond,seg_bonds*2*sizeof(int));
cudaMalloc(&d_angle,seg_angles*3*sizeof(int));
cudaMalloc(&d_dih,seg_dihedrals*4*sizeof(int));
cudaMalloc(&d_impro,seg_impropers*4*sizeof(int));
//cudaMalloc(&d_data,Natom*pd_frame*3*sizeof(float));

cudaMalloc(&d_BondDistance,seg_bonds*seg_frame*sizeof(float));
cudaMalloc(&d_newTheta,seg_angles*seg_frame*sizeof(float));
cudaMalloc(&d_Torsion1,seg_dihedrals*seg_frame*sizeof(float));
cudaMalloc(&d_Torsion2,seg_impropers*seg_frame*sizeof(float));

float *h_BondDistance=(float*) malloc(pd_bonds*pd_frame*sizeof(float));
float *h_newTheta=(float*)malloc(pd_angles*pd_frame*sizeof(float));
float *h_Torsion1=(float*)malloc(pd_dihedrals*pd_frame*sizeof(float));
float *h_Torsion2=(float*)malloc(pd_impropers*pd_frame*sizeof(float));

/////////////////////////Transfer data & call kernell//////////////////////

cudaEventRecord(start, 0);
cudaError_t status = cudaMallocHost((void**)&d_data,Natom*pd_frame*3*sizeof(float));
if (status != cudaSuccess)
   printf("Error allocating pinned host memoryn");
memcpy(d_data,data,Natom*pd_frame*3*sizeof(float));
cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
//////////////////////////////
cudaEventElapsedTime(&time, start, stop);
printf ("Time for the kernel: %f ms\n", time);


cudaStream_t stream1, stream2, stream3, stream4 ;
cudaStreamCreate ( &stream1) ;
cudaStreamCreate ( &stream2) ;
cudaStreamCreate ( &stream3) ;
cudaStreamCreate ( &stream4) ;
int *h_bonds,*h_dihedrals,*h_angles;

//cudaMemcpy(d_bond, bonds, pd_bonds*2*sizeof(int), cudaMemcpyHostToDevice);

HANDLE_ERROR(cudaHostAlloc((void**)&h_bonds, pd_bonds*2*sizeof(int),cudaHostAllocDefault));
for(int v=0;v<pd_bonds*2;v++)    h_bonds[v]=bonds[v]; //check it later if you decrale it at first with out copy is any thing goes wrong with the code

HANDLE_ERROR(cudaHostAlloc((void**)&h_angles, pd_angles*3*sizeof(int),cudaHostAllocDefault));
for(int v=0;v<pd_angles*3;v++)    h_angles[v]=angles[v]; 

HANDLE_ERROR(cudaHostAlloc((void**)&h_dihedrals, pd_dihedrals*4*sizeof(int),cudaHostAllocDefault));
for(int v=0;v<pd_dihedrals*4;v++)    h_dihedrals[v]=dihedrals[v]; 

int* h_impropers;
HANDLE_ERROR(cudaHostAlloc((void**)&h_impropers, pd_impropers*4*sizeof(int),cudaHostAllocDefault));
for(int v=0;v<pd_impropers*4;v++)    h_impropers[v]=impropers[v]; 


cudaEventRecord(start, 0);

while(chunk<FullData){

	if(offset_bonds<pd_bonds) {
		HANDLE_ERROR(cudaMemcpyAsync(d_bond,h_bonds+offset_bonds, seg_bonds*2*sizeof(int), cudaMemcpyHostToDevice,stream1));

		Bdistance<<<dimGridB, dimBlockB,0,stream1>>>(d_bond, d_data,seg_bonds,d_BondDistance);
		cerr << cudaGetErrorString(cudaGetLastError()) << endl;
		err=cudaMemcpyAsync(h_BondDistance+(offset_Distance),d_BondDistance, seg_bonds*seg_frame*sizeof(float), cudaMemcpyDeviceToHost,stream1);
		printf( "%s\n", cudaGetErrorString( err ));

		offset_bonds+=seg_bonds;
		offset_Distance+=seg_bonds*seg_frame;
		chunk+=offset_bonds+offset_Distance; 
		}

	if(offset_angles<pd_angles){
		HANDLE_ERROR(cudaMemcpyAsync(d_angle, h_angles+offset_angles, seg_angles*3*sizeof(int), cudaMemcpyHostToDevice,stream2));
		Angle<<<dimGridA, dimBlockA,0,stream2>>>(d_angle,d_data,seg_angles, d_newTheta);
		cudaMemcpyAsync( h_newTheta+offset_Theta,d_newTheta, seg_angles*seg_frame*sizeof(float), cudaMemcpyDeviceToHost,stream2);
		HANDLE_ERROR(cudaGetLastError());

		offset_angles+=seg_angles;
		offset_Theta+=seg_angles*seg_frame;
		chunk+=offset_angles+offset_Theta;
	 }




 	if(offset_dihedrals<pd_dihedrals){
//for(int h=0;h<pd_angles*3;h++) printf("h_newTheta[%d]=%.10f\n",h,h_newTheta[h]);


		HANDLE_ERROR(cudaMemcpyAsync(d_dih, h_dihedrals+offset_dihedrals, seg_dihedrals*4*sizeof(int), cudaMemcpyHostToDevice,stream3));
		Torsion<<<dimGridP, dimBlockP,0,stream3>>>(d_dih, d_data, seg_dihedrals,d_Torsion1);
		HANDLE_ERROR(cudaGetLastError());
		cudaMemcpyAsync( h_Torsion1+offset_Torsion1,d_Torsion1,seg_dihedrals*seg_frame*sizeof(float), cudaMemcpyDeviceToHost,stream3);

		
		offset_dihedrals+=seg_dihedrals;
		offset_Torsion1+=seg_dihedrals*seg_frame;
		chunk+=offset_dihedrals+offset_Torsion1;

	}



	if(offset_impropers<pd_impropers){
		HANDLE_ERROR(cudaMemcpyAsync(d_impro, h_impropers+offset_impropers, seg_impropers*4*sizeof(int), cudaMemcpyHostToDevice,stream4));
		Torsion<<<dimGridIP, dimBlockIP,0,stream4>>>(d_impro, d_data, seg_impropers,d_Torsion2);
		HANDLE_ERROR(cudaGetLastError());
		cudaMemcpyAsync( h_Torsion2+offset_Torsion2,d_Torsion2,seg_impropers*seg_frame*sizeof(float), cudaMemcpyDeviceToHost,stream4);
		HANDLE_ERROR( cudaStreamSynchronize( stream4 ) );
	       offset_impropers+=seg_impropers;
		offset_Torsion2+=seg_impropers*seg_frame;
		chunk+=offset_impropers+offset_Torsion2;
		}
}
//for(int h=0;h<pd_impropers*4;h++) printf("hTorsion2[%d]=%.10f\n",h,h_Torsion2[h]);
HANDLE_ERROR( cudaStreamSynchronize( stream1 ) );
HANDLE_ERROR( cudaStreamSynchronize( stream2 ) );
HANDLE_ERROR( cudaStreamSynchronize( stream3 ) );
HANDLE_ERROR( cudaStreamSynchronize( stream4 ) );

cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);
//////////////////////////////
cudaEventElapsedTime(&time, start, stop);
printf ("Time for the kernel: %f ms\n", time);
///////***********************
	//free(BondDistance);
  cudaFree(d_bond);
  cudaFree(d_BondDistance);
  cudaFree(d_angle);
  cudaFree(d_newTheta);
  cudaFree(d_dih);
  cudaFree(d_Torsion1);
  cudaFree(d_dih);
  cudaFree(d_Torsion1);


  //cudaFree(d_data);
  return 0;


}









