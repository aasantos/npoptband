#ifndef runfunc_cuh
#define runfunc_cuh
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ls.cuh"
#include "qp.cuh"
#include "rndestim.cuh"
#include "kernelfunc.cuh"
//
//
//
void run_in_cpu(float hcmin,float hcmax,int nhc,float hpmin,float hpmax,int nhp,float xmin,float xmax,int nx)
{
	printf("Start ... \n");
	//
	int nc;
	float *callprice = readArray<float>("callprice.txt",&nc);
	float *callstrike = readArray<float>("callstrike.txt",&nc);
	float *callopenint = readArray<float>("callopenint.txt",&nc);
	printf("nc: %d\n",nc);
	int np;
	float *putprice = readArray<float>("putprice.txt",&np);
	float *putstrike = readArray<float>("putstrike.txt",&np);
	float *putopenint = readArray<float>("putopenint.txt",&np);
	printf("np: %d\n",np);
	float r = 0.02;
	float tau = 1.0/12.0;
	//
	int nu;
	float *strikeunique = readArray<float>("strike.txt",&nu);
	printf("nu: %d\n",nu);
	//
	float *xxx0 = new float[nx];
	float *hc = new float[nhc];
	float *hp = new float[nhp];
	linspace(xmin,xmax,nx,xxx0);
	linspace(hcmin,hcmax,nhc,hc);
	linspace(hpmin,hpmax,nhp,hp);
	//
	int size = nhc*nhp;
	float *cvvector = new float[size];
	float *areavector = new float[size];
	float *entropyvector = new float[size];
	float *varvector = new float[size];
	float *matcrit = new float[size];
	//
	for(int i=0;i<size;i++){
		cvvector[i] = 0.0;
		areavector[i] = 0.0;
		entropyvector[i] = 0.0;
		varvector[i] = 0.0;
		matcrit[i] = 0.0;
	}
	//
	//
	double time_spent = 0.0;
	clock_t begin = clock();
    //
	for(int k=0;k<nu;k++)
	{
		printf("Iteration: %d/%d\n",k+1,nu);
			int iter = -1;
			for(int i=0;i<nhc;i++){
				for(int j=0;j<nhp;j++){
					float cv;
					float area;
					float entropy;
					float variation;
					callputfxex(&cv,&area,&entropy,&variation,
								callprice,callstrike,callopenint,nc,hc[i],
								putprice,putstrike,putopenint,np,hp[j],r,tau,xxx0,nx,strikeunique[k]);
					iter++;
					cvvector[iter] += cv;
					areavector[iter] += area;
					entropyvector[iter] += entropy;
					varvector[iter] += variation;
					matcrit[iter] += cv*variation + (1.0 + fabs(area - 1.0))/entropy;
				}
			}
	}
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("The elapsed time is %f seconds\n", time_spent);
	//
	float hoptim[2];
	minMatrix<float>(hoptim,hc,nhc,hp,nhp,matcrit);
	printf("hc: %.3f; hp: %.3f\n",hoptim[0],hoptim[1]);
	//
	printf("Done ... \n");
	//
	delete[] cvvector;
	delete[] areavector;
	delete[] entropyvector;
	delete[] varvector;
	delete[] matcrit;
	//
	free(callprice);
	free(callstrike);
	free(callopenint);
	free(putprice);
	free(putstrike);
	free(putopenint);
	free(strikeunique);
	//
	delete[] xxx0;
	delete[] hc;
	delete[] hp;
	//
}
//
//
//
void run_in_gpu(float hcmin,float hcmax,int nhc,float hpmin,float hpmax,int nhp,float xmin,float xmax,int nx)
{
	cudaDeviceSetLimit(cudaLimitMallocHeapSize,4194304000L);
	int nc;
	float *callprice = readArray2cuda<float>("callprice.txt",&nc);
	float *callstrike = readArray2cuda<float>("callstrike.txt",&nc);
	float *callopenint = readArray2cuda<float>("callopenint.txt",&nc);
	printf("nc: %d\n",nc);
	int np;
	float *putprice = readArray2cuda<float>("putprice.txt",&np);
	float *putstrike = readArray2cuda<float>("putstrike.txt",&np);
	float *putopenint = readArray2cuda<float>("putopenint.txt",&np);
	printf("np: %d\n",np);
	//
	float r = 0.02;
	float tau = 1.0/12.0;
	//
	int nu;
	float *strikeuniqe = readArray2cuda<float>("strike.txt",&nu);
	printf("nu: %d\n",nu);
	float *xxx0 = linspacecuda<float>(xmin,xmax,nx);
	float *hc = linspacecuda<float>(hcmin,hcmax,nhc);
	float *hp = linspacecuda<float>(hpmin,hpmax,nhp);
	//
	float *areavector;
	cudaMallocManaged(&areavector,nhc*nhp*sizeof(float));
	float *entropyvector;
	cudaMallocManaged(&entropyvector,nhc*nhp*sizeof(float));
	float *varvector;
	cudaMallocManaged(&varvector,nhc*nhp*sizeof(float));
	float *cvvector;
	cudaMallocManaged(&cvvector,nhc*nhp*sizeof(float));
	//
	//  Starting the clock  
	//
	double time_spent = 0.0;
	clock_t begin = clock();
	//
	//
	kernel<<<nhc,nhp>>>(callprice,callstrike,callopenint,nc,hc,
	                      putprice,putstrike,putopenint,np,hp,
	                      xxx0,nx,strikeuniqe,nu,areavector,
	                      entropyvector,varvector,cvvector,r,tau);
	cudaDeviceSynchronize();
	//
	int size = nhc*nhp;
	writeArray<float>(cvvector,"cvvector.txt",size);
	writeArray<float>(areavector,"areavector.txt",size);
	writeArray<float>(varvector,"varvector.txt",size);
	writeArray<float>(entropyvector,"entropyvector.txt",size);
	//
	float *CV = (float*)malloc(size*sizeof(float));
	for(int i=0;i<size;i++) CV[i] = cvvector[i]*varvector[i] + (1.0 + fabs(areavector[i] - 1.0))/entropyvector[i];
	writeArray<float>(CV,"CV.txt",size);
	float hoptim[2];
	minMatrix(hoptim,hc,nhc,hp,nhp,CV);
	printf("hc: %.3f; hp: %.3f\n",hoptim[0],hoptim[1]);
	//
	//
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("The elapsed time is %f seconds\n", time_spent);
	//
	free(CV);
	//
	cudaFree(callprice);
	cudaFree(callstrike);
	cudaFree(callopenint);
	cudaFree(putprice);
	cudaFree(putstrike);
	cudaFree(putopenint);
	cudaFree(strikeuniqe);
	cudaFree(xxx0);
	cudaFree(hc);
	cudaFree(hp);
	cudaFree(areavector);
	cudaFree(entropyvector);
	cudaFree(varvector);
	cudaFree(cvvector);
	printf("Done .... \n");
}

void run_in_gpu_shared(float hcmin,float hcmax,int nhc,float hpmin,float hpmax,int nhp,float xmin,float xmax,int nx)
{
	cudaDeviceSetLimit(cudaLimitMallocHeapSize,4194304000L);
	int nc;
	float *callprice = readArray2cuda<float>("callprice.txt",&nc);
	float *callstrike = readArray2cuda<float>("callstrike.txt",&nc);
	float *callopenint = readArray2cuda<float>("callopenint.txt",&nc);
	printf("nc: %d\n",nc);
	int np;
	float *putprice = readArray2cuda<float>("putprice.txt",&np);
	float *putstrike = readArray2cuda<float>("putstrike.txt",&np);
	float *putopenint = readArray2cuda<float>("putopenint.txt",&np);
	printf("np: %d\n",np);
	//
	float r = 0.02;
	float tau = 1.0/12.0;
	//
	int nu;
	float *strikeuniqe = readArray2cuda<float>("strike.txt",&nu);
	printf("nu: %d\n",nu);
	float *xxx0 = linspacecuda<float>(xmin,xmax,nx);
	float *hc = linspacecuda<float>(hcmin,hcmax,nhc);
	float *hp = linspacecuda<float>(hpmin,hpmax,nhp);
	//
	float *areavector;
	cudaMallocManaged(&areavector,nhc*nhp*sizeof(float));
	float *entropyvector;
	cudaMallocManaged(&entropyvector,nhc*nhp*sizeof(float));
	float *varvector;
	cudaMallocManaged(&varvector,nhc*nhp*sizeof(float));
	float *cvvector;
	cudaMallocManaged(&cvvector,nhc*nhp*sizeof(float));
	//
	//  Starting the clock  
	//
	double time_spent = 0.0;
	clock_t begin = clock();
	//
	//
	kernel_shared<<<dim3(nhc,nhp,1),nx>>>(callprice,callstrike,callopenint,nc,hc,
	                      putprice,putstrike,putopenint,np,hp,
	                      xxx0,nx,strikeuniqe,nu,areavector,
	                      entropyvector,varvector,cvvector,r,tau);
	cudaDeviceSynchronize();
	//
	//
	int size = nhc*nhp;
	writeArray<float>(cvvector,"cvvector.txt",size);
	writeArray<float>(areavector,"areavector.txt",size);
	writeArray<float>(varvector,"varvector.txt",size);
	writeArray<float>(entropyvector,"entropyvector.txt",size);
	//
	float *CV = (float*)malloc(size*sizeof(float));
	for(int i=0;i<size;i++) CV[i] = cvvector[i]*varvector[i] + (1.0 + fabs(areavector[i] - 1.0))/entropyvector[i];
	writeArray<float>(CV,"CV.txt",size);
	float hoptim[2];
	minMatrix(hoptim,hc,nhc,hp,nhp,CV);
	printf("hc: %.3f; hp: %.3f\n",hoptim[0],hoptim[1]);
	//
	//
	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("The elapsed time is %f seconds\n", time_spent);
	//
	free(CV);
	//
	cudaFree(callprice);
	cudaFree(callstrike);
	cudaFree(callopenint);
	cudaFree(putprice);
	cudaFree(putstrike);
	cudaFree(putopenint);
	cudaFree(strikeuniqe);
	cudaFree(xxx0);
	cudaFree(hc);
	cudaFree(hp);
	cudaFree(areavector);
	cudaFree(entropyvector);
	cudaFree(varvector);
	cudaFree(cvvector);
	printf("Done .... \n");
}


#endif