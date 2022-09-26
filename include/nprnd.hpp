//
//  rnd.hpp
//  agosto22
//
//  Created by Antonio Santos on 29/08/2022.
//

#ifndef nprnd_hpp
#define nprnd_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include "ls.hpp"
#include "qp.hpp"
//
//
//
//
class NPRND{
protected:
    int ncall;    // number of observations (call)
    double *callp;  // call prices
    double *calls;  // call strikes
    double *callw;  // call weights
    //
    int nput;    // number of observations (put)
    double *putp;  // put prices
    double *puts; // put strikes
    double *putw; // put weights
    //
    double r;    // interest risk-free rate
    double tau;  // time to maturiry
    //
    //
    double *xRange;
    double *ddf;
    double xRangemin;
    double xRangemax;
    int mRange;     // number of grid elements
    //
    double *strike;   // vector of unique strikes in sample
    int nstrike;     // vector dimension
    //
    //
    // computational statistics
    unsigned long nqpproblems;
    unsigned long niterations;
    //
    double solqp[8];
    //
    //
public:
    //
    NPRND(double *callp,double *calls,double *callw,int ncall,
          double *putp,double *puts,double *putw,int nput,
          double *xRange,int mRange,
          double *strike,int nstrike,
	  double r,double tau)
    {
        //
	//
	//
        this->callp = callp;
        this->calls = calls;
        this->callw = callw;
        this->ncall = ncall;
	//
        this->putp = putp;
        this->puts = puts;
        this->putw = putw;
        this->nput = nput;
       	//
	this->xRange = xRange;
	this->mRange = mRange;
	//
	//
	this->strike = strike;
	this->nstrike = nstrike;
        //
        //
        this->r = r;
        this->tau = tau;
        //        
        //
        this->nqpproblems = 0;
        this->niterations = 0;
    }
    //
    //
    ~NPRND()
    {
	/*
        if(callp) free(callp);
        if(calls) free(calls);
        if(callw) free(callw);
        if(putp) free(putp);
        if(puts) free(puts);
        if(putw) free(putw);
        if(strike) free(strike);
        if(xRange) free(xRange);
        if(ddf) free(ddf);
	*/
    }
    //
    //
    void estim(double x0,double hc,double hp)
    {
        //
        //
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
        for(int i=0;i<this->ncall;i++){
            double t0 = (this->calls[i] - x0);
            double tt = sqrt(this->callw[i]*exp(-0.5*t0*t0/hc/hc));
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
            xy1 += tt*this->callp[i]*t1;
            xy2 += tt*this->callp[i]*t2;
            xy3 += tt*this->callp[i]*t3;
            xy4 += tt*this->callp[i]*t4;
            this->niterations++;
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
            for(int i=0;i<this->nput;i++){
                double t0 = (this->puts[i] - x0);
                double tt = sqrt(this->putw[i]*exp(-0.5*t0*t0/hp/hp));
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
                xyp1 += tt*this->putp[i]*t1;
                xyp2 += tt*this->putp[i]*t2;
                xyp3 += tt*this->putp[i]*t3;
                xyp4 += tt*this->putp[i]*t4;
                this->niterations++;
            }
        //
        //
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
        double tt = exp(-1.0*this->r*this->tau);
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
        QP qp = QP(H,f,8,A,b,8,Aeq,beq,3);
        qp.solve(solqp);
        this->nqpproblems++;
    }
    //
    //
    //
    void estim_no_x(double x0,double x1,double hc,double hp)
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
            for(int i=0;i<this->ncall;i++){
                if(this->calls[i] != x1){
                    double t0 = (this->calls[i] - x0);
                    double tt = sqrt(this->callw[i]*exp(-0.5*t0*t0/hc/hc));
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
                    xy1 += tt*this->callp[i]*t1;
                    xy2 += tt*this->callp[i]*t2;
                    xy3 += tt*this->callp[i]*t3;
                    xy4 += tt*this->callp[i]*t4;
                    this->niterations++;
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
            for(int i=0;i<this->nput;i++){
                if(this->puts[i] != x1){
                    double t0 = (this->puts[i] - x0);
                    double tt = sqrt(this->putw[i]*exp(-0.5*t0*t0/hp/hp));
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
                    xyp1 += tt*this->putp[i]*t1;
                    xyp2 += tt*this->putp[i]*t2;
                    xyp3 += tt*this->putp[i]*t3;
                    xyp4 += tt*this->putp[i]*t4;
                    this->niterations++;
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
            QP qp = QP(H,f,8,A,b,8,Aeq,beq,3);
            qp.solve(solqp);
            this->nqpproblems++;
    }
    //
    //
    //
    void estim_cv_element(double *cv,double *area,double *entropy,double *variation,double xex,double hc,double hp)
    {
        double *ddf = (double*)malloc(this->mRange*sizeof(double));
	    for(int i=0;i<this->mRange;i++){
                this->estim_no_x(this->xRange[i],xex,hc,hp);
                ddf[i] = solqp[2];
                this->niterations++;
            }
        *area = this->areaEstim(this->xRange,ddf,this->mRange);
        *entropy = this->entropyEstim(this->xRange,ddf,this->mRange);
        *variation = this->varEstim(ddf,this->mRange);
	    //
	    free(ddf);
        //
        this->estim_no_x(xex,xex,hc,hp);
        double fcall = solqp[0];
        double fput = solqp[4];
        int nct = 0;
        double cvc = 0.0;
        for(int k=0;k<this->ncall;k++){
            if(this->calls[k] == xex){
                nct++;
                double err = this->callp[k] - fcall;
                cvc += err*err;
                this->niterations++;
                }
            }
        if(nct > 0) cvc /= (double)nct;
        //
        int npt = 0;
        double cvp = 0.0;
        for(int k=0;k<this->nput;k++){
            if(this->puts[k] == xex){
                npt++;
                double err = this->putp[k] - fput;
                cvp += err*err;
                this->niterations++;
            }
        }
        if(npt > 0) cvp /= (double)npt;
        *cv = (cvc + cvp);
    }
    //
    //
    //
    void mat_cv(double *matcv,double *hcv,int nhc,double *hpv,int nhp)
    {
        for(int k=0;k<this->nstrike;k++)
        {
            printf("Iteration: %d/%d\n",k+1,this->nstrike);
            int ii = -1;
             for(int i=0;i<nhc;i++){
                 for(int j=0;j<nhp;j++){
                     ii++;
                     double cv;
                     double area;
                     double entropy;
                     double variation;
                     this->estim_cv_element(&cv,&area,&entropy,&variation,this->strike[k],hcv[i],hpv[j]);
                     matcv[ii] += cv*variation + (1.0 + fabs(area - 1.0))/entropy;
                     this->niterations++;
                 }
            }
        }
    }
    //
    //
    double matCVElement(double hc,double hp,int k)
    {
        double cv;
        double area;
        double entropy;
        double variation;
        this->estim_cv_element(&cv,&area,&entropy,&variation,this->strike[k],hc,hp);
        return cv*variation + (1.0 + fabs(area - 1.0))/entropy;
    }
    
    //
    //
    void optim_bandwidth(double *hoptim,double *hcv,int nhc,double *hpv,int nhp)
    {
        double *matcv = (double*)malloc(nhc*nhp*sizeof(double));
        this->mat_cv(matcv, hcv, nhc, hpv, nhp);
        double sol[2];
        this->minMatrix(sol, hcv, nhc, hpv, nhp, matcv);
        free(matcv);
        hoptim[0] = sol[0];
        hoptim[1] = sol[1];
    }
    //
    //
    long get_nqpprob()
    {
        return this->nqpproblems;
    }
    //
    //
    unsigned long get_niterations()
    {
        return this->niterations;
    }
//
//
private:
//
//
//
void minMatrix(double *res,double *xxx,int nx,double *yyy,int ny,double *mat)
{
    double min = mat[0];
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
    res[0] = xxx[row];
    res[1] = yyy[col];
}
//
//
//
double areaEstim(double *x0,double *yy,int n)
{
     double sum = 0.0;
     for(int i=0;i<n;i++){
        if(i==0){
             sum  += 0.5*(x0[i+1] - x0[i])*yy[i];
             }else if(i==(n-1)){
                sum += 0.5*(x0[i] - x0[i-1])*yy[i];
             }else{
                sum += 0.5*(x0[i+1] - x0[i-1])*yy[i];
            }
     }
    return sum;
}
//
//
//
double entropyEstim(double *x0,double *yy,int n)
{
        double area = this->areaEstim(x0, yy, n);
        double *y0 = (double*)calloc(n,sizeof(double));
        for(int i=0;i<n;i++) y0[i] = yy[i]/area;
        double sum = 0.0;
        for(int i=0;i<n;i++){
            if(i==0){
                    if(y0[i] > 0.0001) sum  += 0.5*(x0[i+1] - x0[i])*y0[i]*log(y0[i]);
                }else if(i==(n-1)){
                    if(y0[i] > 0.0001) sum += 0.5*(x0[i] - x0[i-1])*y0[i]*log(y0[i]);
                }else{
                    if(y0[i] > 0.0001) sum += 0.5*(x0[i+1] - x0[i-1])*y0[i]*log(y0[i]);
                }
        }
        free(y0);
        return -1.0*sum;
}
//
//
//
double varEstim(double *yy,int n)
{
   double sum = 0.0;
   for(int i=1;i<n;i++){
      sum += fabs(yy[i] - yy[i-1]);
   }
   return sum;
}
//
//
};
//


#endif /* rnd_hpp */
