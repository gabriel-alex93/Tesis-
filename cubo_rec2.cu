#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "cuda_fp16.h"


int PrintCu(const int N,int *a);
__global__ void MetNormal(int *a_d,const int N );
__global__ void diagonal(int *a_d ,const int N ,int Nm ,int xi ,int yi,int zi,int n_rec);
__global__ void rec(int *a_d ,const int N ,int Nm ,int xi ,int yi,int zi,int n_rec,const  int  N_CORTE, const int BlockSize){

	int x=blockIdx.x*blockDim.x+threadIdx.x+xi;
	int x0=blockIdx.x*blockDim.x+threadIdx.x;
	int y=blockIdx.y*blockDim.y+threadIdx.y+yi;
	int y0=blockIdx.y*blockDim.y+threadIdx.y;
	int z=blockIdx.z*blockDim.z+threadIdx.z+zi;
	int z0=blockIdx.z*blockDim.z+threadIdx.z;
	int tid1=z*N*N+y*N+x;
	int tid=z0*N*N+y0*N+x0;

	int NmF= (int)ceil((float)((Nm)/2.0f));
	//printf("pspace x y z (%i %i %i)   dominio x y z (%i %i %i)  tid1 dominio %i  tid pspace %i \n",x0, y0, z0, x, y, z,tid1,tid);
 
	//if(Nm>=BlockSize ){	
	if(x < N && y < N && z < N){
		//a_d[tid1]=n_rec;
		a_d[tid1]=1;
		if(x+y >z){
			a_d[tid1]=9;
		}
	}
	//}
	if(Nm<=N_CORTE && tid==0){ //  AQUI NCORTE(pensar en calcular N_corte en base al basico optimo <7>)

		int Nbloques= (Nm+BlockSize-1)/BlockSize;

		dim3 b(BlockSize,BlockSize,BlockSize);
		dim3 g(Nbloques,Nbloques,Nbloques);

		cudaStream_t s1;
		cudaStream_t s2;
		cudaStream_t s3;

		cudaStreamCreateWithFlags(&s1,cudaStreamNonBlocking);
		cudaStreamCreateWithFlags(&s2,cudaStreamNonBlocking);		      
		cudaStreamCreateWithFlags(&s3,cudaStreamNonBlocking);	


		diagonal<<<g,b,0,s1>>>(a_d,N,Nm,xi+Nm,yi,zi,n_rec);
		diagonal<<<g,b,0,s2>>>(a_d,N,Nm,xi,yi+Nm,zi,n_rec);
		diagonal<<<g,b,0,s3>>>(a_d,N,Nm,xi,yi,zi-Nm,n_rec);
		//diagonal kernel nuevo 
		return;
	}
	//busca los nuevos puntos a considerar para el siguiente mapeo
	if( tid==0 && Nm!=0){
		//int  Nbloques=(((int) ceil((float)(Nm)/2.0f))+BlockSize-1)/BlockSize;
		int Nbloques=(NmF+BlockSize-1)/BlockSize;
		//printf("nbloques= %i  %i \n", Nbloques, BlockSize);
		dim3 b(BlockSize,BlockSize,BlockSize);
		dim3 g(Nbloques,Nbloques,Nbloques);

		n_rec=n_rec+1;		

 		cudaStream_t s1;
		cudaStream_t s2;
		cudaStream_t s3;

		cudaStreamCreateWithFlags(&s1,cudaStreamNonBlocking);
		cudaStreamCreateWithFlags(&s2,cudaStreamNonBlocking);		      
		cudaStreamCreateWithFlags(&s3,cudaStreamNonBlocking);		
		//
		//int NmF= (int)ceil((float)((Nm)/2.0f));
		//printf("kernel Nm = %i    NmF= %i \n",Nm, NmF);
	
		rec<<<g,b,0,s1>>>(a_d,N,NmF,xi+Nm,yi,zi+NmF,n_rec,N_CORTE,BlockSize);
		rec<<<g,b,0,s2>>>(a_d,N,NmF,xi,yi+Nm,zi+NmF,n_rec,N_CORTE,BlockSize);
		rec<<<g,b,0,s3>>>(a_d,N,NmF,xi,yi,zi-NmF,n_rec,N_CORTE,BlockSize);
	}
}
__global__ void diagonal(int *a_d ,const int N ,int Nm ,int xi ,int yi,int zi,int n_rec){
	int x= blockIdx.x*blockDim.x+threadIdx.x+xi;
	int y= blockIdx.y*blockDim.y+threadIdx.y+yi;
	int z=  blockIdx.z*blockDim.z+threadIdx.z+zi;
	int tid1= z*N*N+y*N+x;
	if(x < N && y < N && z < N){
		if(x+y<= z){
			a_d[tid1]=1;
		}	
	}

}
__global__ void MetNormal(int *a_d,const int N){
	int x= blockIdx.x*blockDim.x+threadIdx.x;
	int y= blockIdx.y*blockDim.y+threadIdx.y;
	int z= blockIdx.z*blockDim.z+threadIdx.z;
	int ind=z*N*N+y*N+x;
	if( x+y<=z){
		a_d[ind]=1;
	}

}
int PrintCu(const int N, int *a){

	for(int k=0;k<N;k++){
		printf("z= %i \n",k);
			for (int i=0;i<N;i++){
				for(int j=0;j<N;j++){
					if(a[N*N*k+i*N+j] != 9){
					printf("%i ",a[N*N*k+i*N+j]);
			   		 }
				    	else{
						printf("* ");
					    }
					}
				printf("\n");
			 }
			 printf("\n");	
		}
	//printf("Ok\n");
	return 0;
}

int Ver_Resultado(const int N, int *a){
	for(int i=0;i<N;i++){
    	for (int j=0;j<N;j++){
        	for(int k=0;k<N;k++){
         		 if(  !((i+j<=k && a[N*N*k+i*N+j]==1) || (i+j>k && a[N*N*k+i*N+j]==9) ) ){
       	  			printf("Error en Matriz \n" );
            	  	exit(1);
              		
            	  }
    		}
    	}
	}
	    //printf("Matriz correcta \n");
        return (0);
}
int main(int argc ,char **argv){
	
	if (argc !=6){
		fprintf(stderr,"error, ejecutar  programa como ./prog N met rep ncorte  BlockSize\n");
		exit(EXIT_FAILURE);
	}
	unsigned long N=atoi(argv[1]);
	unsigned long met=atoi(argv[2]);
	unsigned long rep=atoi(argv[3]);
	unsigned long ncort=atoi(argv[4]);
	//int  nt=atoi(argv[5]);//numero de threads por bloque(eliminar )
	int BSize=atoi(argv[5]);
	int *a,*a_d, xi=0,yi=0,zi=(int) ceil((float)(N/2.0f));
	//double *datos;

	//printf("malloc ...");
	fflush(stdout);
	a=(int*)malloc(sizeof(int)*N*N*N);
	//datos=(double*)malloc(sizeof(double)*rep);

	//printf("ok ...\ncuda malloc...");
	fflush(stdout);
	cudaMalloc((void ** ) &a_d,N*N*N*sizeof(int));
	//printf("ok ...\n");
	fflush(stdout);
	
	dim3 Bloque(BSize,BSize,BSize);//un  bloquede nt 

	float NB=(float)N/(float)(2*BSize);
	int B=(int) ceil(NB);
	dim3 Grid(B,B,B);//bgrid  de B*b*b bloque
	dim3 GridBruto((N+BSize-1)/BSize,(N+BSize-1)/BSize,(N+BSize-1)/BSize);
	//printf("inicializando con N= %i ...",N);
	fflush(stdout);
	for(int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			for(int k=0;k<N;k++){
				a[N*N*k+i*N+j]=9;
			}
		}
	}
	//printf(" ok..\n");
	fflush(stdout);
	int n_rec=0;	
	double t1=omp_get_wtime();
	cudaMemcpy(a_d,a,N*N*N*sizeof(int),cudaMemcpyHostToDevice);
	//printf("calculo GPU...\n");
	fflush(stdout);
	double t2;
	double t3;
	if(ncort >= BSize && (met==1 || met==2 )){
		if(met==1){// aqui se supone que viene un while o for para las iteraciones 
			//printf("Metodo recursivo......\n"); fflush(stdout);
			for(int i=0;i<150;i++){
				rec<<<Grid,Bloque>>>(a_d,N,(int) ceil((float)(N)/2.0f),xi,yi,zi,n_rec,ncort,BSize);
 				cudaDeviceSynchronize();	
 				cudaError_t error = cudaGetLastError();
  				if(error != cudaSuccess)
  				{
    // print the CUDA error message and exit
    			printf("CUDA error: %s\n", cudaGetErrorString(error));
    			exit(-1);
  }
			}
			t2=omp_get_wtime();
			for(int i=0;i<rep;i++){
				rec<<<Grid,Bloque>>>(a_d,N,(int) ceil((float)(N)/2.0f),xi,yi,zi,n_rec,ncort,BSize);
 				cudaDeviceSynchronize();	
 				cudaError_t error = cudaGetLastError();
  				if(error != cudaSuccess)
  				{
  				  // print the CUDA error message and exit
   				 printf("CUDA error: %s\n", cudaGetErrorString(error));
    			exit(-1);
 				 }
			}
			t3=omp_get_wtime();
			//printf("ok\n"); 
			fflush(stdout);
		}
		if(met==2){
			//printf("Metodo bruto...\n");
			for(int i =0;i<150;i++){
				MetNormal<<<GridBruto,Bloque>>>(a_d,N);	
				cudaDeviceSynchronize();
				}
			t2=omp_get_wtime();
			for(int i =0;i<rep;i++){
				MetNormal<<<GridBruto,Bloque>>>(a_d,N);	
				cudaDeviceSynchronize();
				}
			t3=omp_get_wtime();
			fflush(stdout);
			}
		
	}
	else{
		printf("Error, N de corte menor a tamaño de bloque  o metodo invalido\n");
		return;
	}
	//aqui calculo el promedio de los tiempos
	double media=(t3-t2)/rep;
	//printf("Tiempo promedio con %i iteraciones: %f \n",rep,media);
	fflush(stdout);

	cudaDeviceSynchronize();
	//printf("ok..\n");
	cudaMemcpy(a,a_d,N*N*N*sizeof(int),cudaMemcpyDeviceToHost);
	
	double t4=omp_get_wtime();
	/*printf("calculo cpu...");
	fflush(stdout);
	double t5=omp_get_wtime();
	printf("ok..\n");
	printf("verificando...\n");
    */
   	 
	//if(N < 128){
		//PrintCu(N,a);//imprime cubo
   	//}
		Ver_Resultado(N,a);
	//printf("grid : %i %i %i,Bloque:  %i  %i  %i \n",Bloque.x,Bloque.y,Bloque.z ,Grid.x,Grid.y,Grid.z);
	//printf("gridBruto : %i %i %i,Bloque:  %i  %i  %i \n",GridBruto.x,GridBruto.y,GridBruto.z ,Grid.x,Grid.y,Grid.z);
	//printf("tiempo copy a gpu : %f\ntiempo kernel: %f\ntiempo copy to host: %f tiempo total: %f\n",t2-t1,media,t4-t3,t4-t1);
	//printf("tiempo cpu %f\n",t5-t4);
	printf("%f\n",1000*media);
	return 0;
	}

