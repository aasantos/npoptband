#ifndef runfunc_hpp
#define runfunc_hpp



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include "ls.hpp"
#include "qp.hpp"
#include "nprnd.hpp"


template <typename T>
T* readArray(const char* file,int *n);
//
//
void minMatrix(double *res,double *xxx,int nx,double *yyy,int ny,double *mat);
double* linspace(double x0,double x1,int n);
double areaEstim(double *x0,double *yy,int n);
double entropyEstim(double *x0,double *yy,int n);
double varEstim(double *yy,int n);
//
// needs to define the grids for hc (hcmin,hcmax,nhc) and hp (hpmin,hpmax,nhp)
// and the grid for the support of the density (xmin,xmax,nx)
// define the risk-free interest rate (r) and time to maturity 
void run_in_cpu(double hcmin,double hcmax,int nhc,double hpmin,double hpmax,int nhp,double xmin,double xmax,int nx,double r,double tau);


#endif
