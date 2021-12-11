#ifndef kernelfunc_cuh
#define kernelfunc_cuh
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ls.cuh"
#include "qp.cuh"
#include "rndestim.cuh"
//
//
__global__ void kernel(float *y,float *x,float *w,int nc,float *hc,
						 float *z,float *zx,float *zw,int np,float *hp,
						 float *xxx0,int nx,float *strikeunique,int nu,
						 float *areavector,float *entropyvector,float *varvector,float *cvvector,
						 float r,float tau)
{
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	//
	float *areavec = (float*)malloc(nu*sizeof(float));
	float *entropyvec = (float*)malloc(nu*sizeof(float));
	float *varvec = (float*)malloc(nu*sizeof(float));
	float *cvvec = (float*)malloc(nu*sizeof(float));
	for(int i=0;i<nu;i++)
	{
		float *rnd_temp = (float*)malloc(nx*sizeof(float));
		for(int j=0;j<nx;j++)
		{
			float sol[8];
			npcallputoptim_no_x1(sol,y,x,w,nc,hc[blockIdx.x],z,zx,zw,np,hp[threadIdx.x],
			xxx0[j],strikeunique[i],r,tau);
			rnd_temp[j] = sol[2];
		}
		float tt = area_estim(xxx0,rnd_temp,nx);
		float tt1 = entropy_estim(xxx0,rnd_temp,nx);
		float tt2 = variation_estim(rnd_temp,nx);
		areavec[i] = tt;
		entropyvec[i] = tt1;
		varvec[i] = tt2;
		cvvec[i] = cv_estim(y,x,w,nc,hc[blockIdx.x],z,zx,zw,np,hp[threadIdx.x],strikeunique[i],r,tau);
		free(rnd_temp);
	}
	//
	float aa = mean(areavec,nu);
	float aa1 = mean(entropyvec,nu);
	float aa2 = mean(varvec,nu);
	float aa3 = mean(cvvec,nu);
	free(areavec);
	free(entropyvec);
	free(varvec);
	free(cvvec);
	areavector[idx] = aa;
	entropyvector[idx] = aa1;
	varvector[idx] = aa2;
	cvvector[idx] = aa3;
}
//////////////////////////////////////
__global__ void kernel_shared(float *y,float *x,float *w,int nc,float *hc,
						 float *z,float *zx,float *zw,int np,float *hp,
						 float *xxx0,int nx,float *strikeunique,int nu,
						 float *areavector,float *entropyvector,float *varvector,
						 float *cvvector,float r,float tau)
{
	int idx = gridDim.x*blockIdx.x + blockIdx.y;
	//
	float *areavec = (float*)malloc(nu*sizeof(float));
	float *entropyvec = (float*)malloc(nu*sizeof(float));
	float *varvec = (float*)malloc(nu*sizeof(float));
	float *cvvec = (float*)malloc(nu*sizeof(float));
	//
	__shared__ float rnd_temp[128];
	float sol[8];
	//
	for(int i=0;i<nu;i++)
	{
			//
			npcallputoptim_no_x1(sol,y,x,w,nc,hc[blockIdx.x],z,zx,zw,np,hp[blockIdx.y],
			xxx0[threadIdx.x],strikeunique[i],r,tau);
			rnd_temp[threadIdx.x]  = sol[2];
			__syncthreads();
			//
			float tt = area_estim(xxx0,rnd_temp,nx);
			float tt1 = entropy_estim(xxx0,rnd_temp,nx);
			float tt2 = variation_estim(rnd_temp,nx);
			float tt3 = cv_estim(y,x,w,nc,hc[blockIdx.x],z,zx,zw,np,hp[blockIdx.y],strikeunique[i],r,tau);
			//
			areavec[i] = tt;
			entropyvec[i] = tt1;
			varvec[i] = tt2;
			cvvec[i] = tt3;
			//
	}
	//
	float aa = mean(areavec,nu);
	float aa1 = mean(entropyvec,nu);
	float aa2 = mean(varvec,nu);
	float aa3 = mean(cvvec,nu);
	//
	free(areavec);
	free(entropyvec);
	free(varvec);
	free(cvvec);
	//
	areavector[idx] = aa;
	entropyvector[idx] = aa1;
	varvector[idx] = aa2;
	cvvector[idx] = aa3;
}



#endif
