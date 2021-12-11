#ifndef vector_cuh
#define vector_cuh

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>


//////////////////////////////////////////////////////////////////////
__host__ __device__ void linspace(float x0,float x1,int n,float *xx)
{
	float step = (x1 - x0)/(float)(n - 1);
	xx[0] = x0;
	for(int i=1;i<n;i++) xx[i] = xx[i-1] + step;
}
//
//
__host__ __device__ float* linspace(float x0,float x1,int n)
{
	float *xx = (float*)malloc(n*sizeof(float));
	float step = (x1 - x0)/(float)(n - 1);
	xx[0] = x0;
	for(int i=1;i<n;i++) xx[i] = xx[i-1] + step;
	return xx;
}
//
//
template <typename T>
__host__ T* linspacecuda(T x0,T x1,int n)
{
	T* xx;
	cudaMallocManaged(&xx,n*sizeof(T));
	float step = (x1 - x0)/(float)(n - 1);
	xx[0] = x0;
	for(int i=1;i<n;i++) xx[i] = xx[i-1] + step;
	return xx;
}
///////////////////////////////////////////////////////////////////////

__host__ __device__ float mean(float *x,int n)
{
	float sum = 0.0;
	for(int i=0;i<n;i++) sum += x[i];
	return sum/(float)n;
}
//
//
template <typename T>
__host__ __device__ T sum(T*x,int n){
    T sum = 0.0;
    for(int i=0;i<n;i++) sum += x[i];
    return sum;
}
//
//
template <typename  T>
void minMatrix(T* res,T *x,int nx,T *y,int ny,T *mat)
{
    T min = mat[0];
    int row = 0;
    int col = 0;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            if(mat[i*nx + j] < min){
            	min = mat[i*nx + j];
                row = i;
                col = j;
            }
        }
    }
    res[0] = x[row];
    res[1] = y[col];
}
//
//
template <typename T>
	void swap(T *a, T *b)
{
    T temp = *a;
    *a = *b;
    *b = temp;
}

template <typename T>
void bubbleSort(T array[], int n)
{
    for (int i = 0; i < n-1; i++)
        for (int j = 0; j < n-i-1; j++) if (array[j] > array[j+1])
            swap<T>(&array[j], &array[j+1]);
}
//
template <typename T>
void sort(T *x,int n)
{
    bubbleSort<T>(x,n);
}
//
//
template <typename T>
T* unique(T *x,int n,int *nu)
{
	int nunique = 1;
	sort<T>(x,n);
	for(int i=1;i<n;i++) if(x[i] > x[i-1]) nunique++;
	*nu = nunique;
	T *xout;
	cudaMallocManaged(&xout,nunique*sizeof(T));
	xout[0] = x[0];
	int iter = 0;
	for(int i=1;i<n;i++){
		if(x[i] > x[i-1]){
			iter++;
			xout[iter] = x[i];
		}
	}
	return xout;
}


#endif
