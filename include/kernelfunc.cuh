#ifndef kernelfunc_cuh
#define kernelfunc_cuh

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "nprnd.cuh"
//
/*
*out output vector with dimension (nhp x nhc x nu)


//
*hc *hp *strike input vectors defining the grid and parallel components

// sample data
*callprice *callstrike *callopenint nc  -  call data
*putprice *putstrike *putopenint np     -  put data

// density support components
*xRange mRange 

// option contract parameters
tau -  time to maturity
r   -  risk-free interest rate

*/
//
__global__ void kernel(double *out,
		       double *hc,int nhc,double *hp,int nhp,double *strike,int nu,
                       double *callprice,double *callstrike,double *callopenint,int nc,
                       double *putprice,double *putstrike,double *putopenint,int np,
                       double *xRange,int mRange,
                       double tau,double r)
{
        int ncc = gridDim.x;
        int npp = gridDim.y;
        int xi = blockIdx.x;
        int yi = blockIdx.y;
        int zi = threadIdx.x;
        unsigned int size = nhp*nhc*nu;
        int idx = zi*(ncc*npp - 1) + (ncc*xi + yi) + zi;
        //
        if(idx < size){
        NPRND nprnd = NPRND(callprice,callstrike,callopenint,nc,
                             putprice,putstrike,putopenint,np,
                             r,tau,
                             xRange,mRange,
                             strike,nu);
        out[idx] = nprnd.matCVElement(hc[xi],hp[yi],zi);
        }
}



#endif
