#ifndef rndestim_cuh
#define rndestim_cuh

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "ls.cuh"
#include "qp.cuh"
//
//
//
__host__  __device__ float  area_estim(float *x0,float *y,int n)
{
	float sum = 0.0;
	for(int i=0;i<n;i++){
		if(i==0){
                sum  += 0.5*(x0[i+1] - x0[i])*y[i];
            }else if(i==(n-1)){
                sum += 0.5*(x0[i] - x0[i-1])*y[i];
            }else{
                sum += 0.5*(x0[i+1] - x0[i-1])*y[i];
            }
	}
	return sum;
}
//
//
//
template <typename T>
__host__  __device__ T entropy_estim(T *x0,T *y,int n)
{
    float area = area_estim(x0, y, n);
    float *y0 = (float*)malloc(n*sizeof(float));
    for(int i=0;i<n;i++) y0[i] = y[i]/area;
    float sum = 0.0;
    for(int i=0;i<n;i++){
        if(i==0){
                sum  += 0.5*(x0[i+1] - x0[i])*y0[i]*log(y0[i] + 0.000001);
            }else if(i==(n-1)){
                sum += 0.5*(x0[i] - x0[i-1])*y0[i]*log(y0[i] + 0.000001);
            }else{
                sum += 0.5*(x0[i+1] - x0[i-1])*y0[i]*log(y0[i] + 0.000001);
            }
    }
    free(y0);
    return -1.0*sum;
}
//
//
template <typename T>
__host__ __device__ T variation_estim(T *y,int n)
{
    float sum = 0.0;
    for(int i=1;i<n;i++) sum += fabs(y[i] - y[i-1]);
    return sum;
}
//
//
//
__host__ __device__ void npcallputoptim(float *res,
										float *y,float *x,float *w,int nc,float hc,
										float *z,float *zx,float *zw,int np,float hp,
										float x0,float r,float tau)
{
		double x11 = 0.0; double x12 = 0.0; double x13 = 0.0; double x14 = 0.0;
		double x22 = 0.0; double x23 = 0.0; double x24 = 0.0;
		double x33 = 0.0; double x34 = 0.0;
		double x44 = 0.0;
		//
		double xy1 = 0.0;
		double xy2 = 0.0;
		double xy3 = 0.0;
		double xy4 = 0.0;
		//
		for(int i=0;i<nc;i++){
			double t0 = (x[i] - x0);
			double tt = sqrt(w[i]*exp(-0.5*t0*t0/(hc*hc)));
			double t1 = tt;
			double t2 = t1*t0;
			double t3 = 0.5*t2*t0;
			double t4 = (2.0/6.0)*t3*t0;
			x11 += t1*t1;
			x12 += t1*t2;
			x13 += t1*t3;
			x14 += t1*t4;
			x22 += t2*t2;
			x23 += t2*t3;
			x24 += t2*t4;
			x33 += t3*t3;
			x34 += t3*t4;
			x44 += t4*t4;
			xy1 += tt*y[i]*t1;
			xy2 += tt*y[i]*t2;
			xy3 += tt*y[i]*t3;
			xy4 += tt*y[i]*t4;
		}
		double xp11 = 0.0; double xp12 = 0.0; double xp13 = 0.0; double xp14 = 0.0;
		double xp22 = 0.0; double xp23 = 0.0; double xp24 = 0.0;
		double xp33 = 0.0; double xp34 = 0.0;
		double xp44 = 0.0;
		//
		double xyp1 = 0.0;
		double xyp2 = 0.0;
		double xyp3 = 0.0;
		double xyp4 = 0.0;
		//
			for(int i=0;i<np;i++){
				double t0 = (zx[i] - x0);
				double tt = sqrt(zw[i]*exp(-0.5*t0*t0/(hp*hp)));
				double t1 = tt;
				double t2 = t1*t0;
				double t3 = 0.5*t2*t0;
				double t4 = (2.0/6.0)*t3*t0;
				xp11 += t1*t1;
				xp12 += t1*t2;
				xp13 += t1*t3;
				xp14 += t1*t4;
				xp22 += t2*t2;
				xp23 += t2*t3;
				xp24 += t2*t4;
				xp33 += t3*t3;
				xp34 += t3*t4;
				xp44 += t4*t4;
				//
				xyp1 += tt*z[i]*t1;
				xyp2 += tt*z[i]*t2;
				xyp3 += tt*z[i]*t3;
				xyp4 += tt*z[i]*t4;
			}

		double H[64];
		H[0] = x11;  H[1]  = x12; H[2]  = x13; H[3] =  x14; H[4] =  0.0;  H[5] =  0.0;  H[6]  = 0.0;  H[7] =  0.0;
		H[8] = x12;  H[9]  = x22; H[10] = x23; H[11] = x24; H[12] = 0.0;  H[13] = 0.0;  H[14] = 0.0;  H[15] = 0.0;
		H[16] = x13; H[17] = x23; H[18] = x33; H[19] = x34; H[20] = 0.0;  H[21] = 0.0;  H[22] = 0.0;  H[23] = 0.0;
		H[24] = x14; H[25] = x24; H[26] = x34; H[27] = x44; H[28] = 0.0;  H[29] = 0.0;  H[30] = 0.0;  H[31] = 0.0;
		H[32] = 0.0; H[33] = 0.0; H[34] = 0.0; H[35] = 0.0; H[36] = xp11; H[37] = xp12; H[38] = xp13; H[39] = xp14;
		H[40] = 0.0; H[41] = 0.0; H[42] = 0.0; H[43] = 0.0; H[44] = xp12; H[45] = xp22; H[46] = xp23; H[47] = xp24;
		H[48] = 0.0; H[49] = 0.0; H[50] = 0.0; H[51] = 0.0; H[52] = xp13; H[53] = xp23; H[54] = xp33; H[55] = xp34;
		H[56] = 0.0; H[57] = 0.0; H[58] = 0.0; H[59] = 0.0; H[60] = xp14; H[61] = xp24; H[62] = xp34; H[63] = xp44;
		//
		//
		double f[8];
		f[0] = -1.0*xy1;
		f[1] = -1.0*xy2;
		f[2] = -1.0*xy3;
		f[3] = -1.0*xy4;
		f[4] = -1.0*xyp1;
		f[5] = -1.0*xyp2;
		f[6] = -1.0*xyp3;
		f[7] = -1.0*xyp4;
		//
		double tt = exp(-1.0*r*tau);
        double A[64] = {1,0,0,0,0,0,0,0,
                        0,1,-1,0,0,0,0,0,
                        0,0,0,1,0,0,0,0,
                        0,0,0,0,0,0,0,0,
                        0,0,0,0,1,0,0,0,
                        0,0,0,0,0,1,1,0,
                        0,0,0,0,0,0,0,1,
                        0,0,0,0,0,0,0,0
                        };
        //
        double b[8] = {0,-1.0*tt,0,0,0,0,-1.0*tt,0};
        //
        double Aeq[24] =  {0,0,0,
                           0,0,-1,
                           0,-1,0,
                           -1,0,0,
                            0,0,0,
                            0,0,1,
                            0,1,0,
                            1,0,0};
        //
        double beq[3] = {0,0,tt};
		//
		double sol[8];
		QP *qp = new QP(H,f,8,A,b,8,Aeq,beq,3);
		qp->solve(sol);
		for(int i=0;i<8;i++) res[i] = (float)sol[i];
		delete qp;
}
//
//
//
__host__ __device__ void npcallputoptim_no_x1(float *res,
										float *y,float *x,float *w,int nc,float hc,
										float *z,float *zx,float *zw,int np,float hp,
										float x0,float x1,float r,float tau)
{
		double x11 = 0.0; double x12 = 0.0; double x13 = 0.0; double x14 = 0.0;
		double x22 = 0.0; double x23 = 0.0; double x24 = 0.0;
		double x33 = 0.0; double x34 = 0.0;
		double x44 = 0.0;
		//
		double xy1 = 0.0;
		double xy2 = 0.0;
		double xy3 = 0.0;
		double xy4 = 0.0;
		//
		for(int i=0;i<nc;i++){
			if(x[i] != x1){
				double t0 = (x[i] - x0);
				double tt = sqrt(w[i]*exp(-0.5*t0*t0/(hc*hc)));
				double t1 = tt;
				double t2 = t1*t0;
				double t3 = 0.5*t2*t0;
				double t4 = (2.0/6.0)*t3*t0;
				x11 += t1*t1;
				x12 += t1*t2;
				x13 += t1*t3;
				x14 += t1*t4;
				x22 += t2*t2;
				x23 += t2*t3;
				x24 += t2*t4;
				x33 += t3*t3;
				x34 += t3*t4;
				x44 += t4*t4;
				xy1 += tt*y[i]*t1;
				xy2 += tt*y[i]*t2;
				xy3 += tt*y[i]*t3;
				xy4 += tt*y[i]*t4;
			}

		}
		double xp11 = 0.0; double xp12 = 0.0; double xp13 = 0.0; double xp14 = 0.0;
		double xp22 = 0.0; double xp23 = 0.0; double xp24 = 0.0;
		double xp33 = 0.0; double xp34 = 0.0;
		double xp44 = 0.0;
		//
		double xyp1 = 0.0;
		double xyp2 = 0.0;
		double xyp3 = 0.0;
		double xyp4 = 0.0;
		//
			for(int i=0;i<np;i++){
				if(zx[i] != x1){
					double t0 = (zx[i] - x0);
					double tt = sqrt(zw[i]*exp(-0.5*t0*t0/(hp*hp)));
					double t1 = tt;
					double t2 = t1*t0;
					double t3 = 0.5*t2*t0;
					double t4 = (2.0/6.0)*t3*t0;
					xp11 += t1*t1;
					xp12 += t1*t2;
					xp13 += t1*t3;
					xp14 += t1*t4;
					xp22 += t2*t2;
					xp23 += t2*t3;
					xp24 += t2*t4;
					xp33 += t3*t3;
					xp34 += t3*t4;
					xp44 += t4*t4;
					//
					xyp1 += tt*z[i]*t1;
					xyp2 += tt*z[i]*t2;
					xyp3 += tt*z[i]*t3;
					xyp4 += tt*z[i]*t4;
				}

			}

		double H[64];
		H[0] = x11;  H[1]  = x12; H[2]  = x13; H[3] =  x14; H[4] =  0.0;  H[5] =  0.0;  H[6]  = 0.0;  H[7] =  0.0;
		H[8] = x12;  H[9]  = x22; H[10] = x23; H[11] = x24; H[12] = 0.0;  H[13] = 0.0;  H[14] = 0.0;  H[15] = 0.0;
		H[16] = x13; H[17] = x23; H[18] = x33; H[19] = x34; H[20] = 0.0;  H[21] = 0.0;  H[22] = 0.0;  H[23] = 0.0;
		H[24] = x14; H[25] = x24; H[26] = x34; H[27] = x44; H[28] = 0.0;  H[29] = 0.0;  H[30] = 0.0;  H[31] = 0.0;
		H[32] = 0.0; H[33] = 0.0; H[34] = 0.0; H[35] = 0.0; H[36] = xp11; H[37] = xp12; H[38] = xp13; H[39] = xp14;
		H[40] = 0.0; H[41] = 0.0; H[42] = 0.0; H[43] = 0.0; H[44] = xp12; H[45] = xp22; H[46] = xp23; H[47] = xp24;
		H[48] = 0.0; H[49] = 0.0; H[50] = 0.0; H[51] = 0.0; H[52] = xp13; H[53] = xp23; H[54] = xp33; H[55] = xp34;
		H[56] = 0.0; H[57] = 0.0; H[58] = 0.0; H[59] = 0.0; H[60] = xp14; H[61] = xp24; H[62] = xp34; H[63] = xp44;
		//
		//
		double f[8];
		f[0] = -1.0*xy1;
		f[1] = -1.0*xy2;
		f[2] = -1.0*xy3;
		f[3] = -1.0*xy4;
		f[4] = -1.0*xyp1;
		f[5] = -1.0*xyp2;
		f[6] = -1.0*xyp3;
		f[7] = -1.0*xyp4;
		//
		double tt = exp(-1.0*r*tau);
        double A[64] = {1,0,0,0,0,0,0,0,
                        0,1,-1,0,0,0,0,0,
                        0,0,0,1,0,0,0,0,
                        0,0,0,0,0,0,0,0,
                        0,0,0,0,1,0,0,0,
                        0,0,0,0,0,1,1,0,
                        0,0,0,0,0,0,0,1,
                        0,0,0,0,0,0,0,0
                        };
        //
        double b[8] = {0,-1.0*tt,0,0,0,0,-1.0*tt,0};
        //
        double Aeq[24] =  {0,0,0,
                           0,0,-1,
                           0,-1,0,
                           -1,0,0,
                            0,0,0,
                            0,0,1,
                            0,1,0,
                            1,0,0};
        //
        double beq[3] = {0,0,tt};
		//
		double sol[8];
		QP *qp = new QP(H,f,8,A,b,8,Aeq,beq,3);
		qp->solve(sol);
		for(int i=0;i<8;i++) res[i] = (float)sol[i];
		delete qp;
}
//
//
__host__ __device__ void callputf(float *cf,float *pf,float *ddf,float *area,float *entropy,float *variation,
								  float *y,float *x,float *w,int nc,float hc,
								  float *z,float *zx,float *zw,int np,float hp,
								  float r,float tau,float *xx,int m)
{
	for(int i=0;i<m;i++){
		float sol[8];
		npcallputoptim(sol,y,x,w,nc,hc,z,zx,zw,np,hp,xx[i],r,tau);
		cf[i] = sol[0];
		pf[i] = sol[4];
		ddf[i] = sol[2];
	}
	*area = area_estim(xx,ddf,m);
	*entropy = entropy_estim<float>(xx,ddf,m);
	*variation = variation_estim<float>(ddf,m);
}

__host__ __device__ void callputfxex(float *cv,float *area,float *entropy,float *variation,
								     float *y,float *x,float *w,int nc,float hc,
								     float *z,float *zx,float *zw,int np,float hp,
								     float r,float tau,float *xx,int m,float xex)
{
	float *ddf = new float[m];
	for(int i=0;i<m;i++){
		float sol[8];
		npcallputoptim_no_x1(sol,y,x,w,nc,hc,z,zx,zw,np,hp,xx[i],xex,r,tau);
		ddf[i] = sol[2];
	}
	*area = area_estim(xx,ddf,m);
	*entropy = entropy_estim<float>(xx,ddf,m);
	*variation = variation_estim<float>(ddf,m);
	//
	float sol[8];
	npcallputoptim_no_x1(sol,y,x,w,nc,hc,z,zx,zw,np,hp,xex,xex,r,tau);
	float fcall = sol[0];
	float fput = sol[4];
	int nct = 0;
	float cvc = 0.0;
	for(int k=0;k<nc;k++){
		if(x[k] == xex){
			nct++;
			float err = y[k] - fcall;
			cvc += err*err;
		}
	}
   if(nct > 0) cvc /= (float)nct;
	//
	int npt = 0;
	float cvp = 0.0;
	for(int k=0;k<np;k++){
		if(zx[k] == xex){
			npt++;
			float err = z[k] - fput;
			cvp += err*err;
		}
	}
	if(npt > 0) cvp /= (float)npt;
	*cv = (cvc + cvp);
	delete[] ddf;
}
//
//
//
__device__ float cv_estim(float *y,float *x,float *w,int nc,float hc,
                	      float *z,float *zx,float *zw,int np,float hp,
                          float x0,float r,float tau)
{
    float sol[8];
    npcallputoptim_no_x1(sol,y,x,w,nc,hc,z,zx,zw,np,hp,x0,x0,r,tau);
    float fcall = sol[0];
    float fput = sol[4];
    int nct = 0;
    float cvc = 0.0;
    for(int m=0;m<nc;m++){
        if(x[m] == x0){
            nct++;
            float err = y[m] - fcall;
            cvc += err*err;
        }
    }
   if(nct > 0) cvc /= (float)nct;
    //
    int npt = 0;
    float cvp = 0.0;
    for(int m=0;m<np;m++){
        if(zx[m] == x0){
            npt++;
            float err = z[m] - fput;
            cvp += err*err;
        }
    }
    if(npt > 0) cvp /= (float)npt;
    return (cvc + cvp);
}
//
//


#endif
