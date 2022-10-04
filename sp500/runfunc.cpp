#include "runfunc.hpp"
//
//
template <typename T>
T* readArray(const char* file,int *n)
{
        FILE *fp;
        fp = fopen(file, "r");
        char buff[64];
        int nrow = -1;
        while (!feof(fp)) {
            int tt = fscanf(fp, "%s",buff);
            nrow++;
        }
        *n = nrow;
        rewind(fp);
        T *result = (T*)malloc(nrow*sizeof(T));
        for(int i=0;i<nrow;i++){
            int tt = fscanf(fp, "%s",buff);
            result[i] = (T)atof(buff);
        }
        fclose(fp);
        return result;
}
//
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
//
double* linspace(double x0,double x1,int n)
{
    double *xxx = (double*)malloc(n*sizeof(double));
    double step = (x1 - x0)/(double)(n - 1);
    xxx[0] = x0;
    for(int i=1;i<n;i++){
        xxx[i] = xxx[i-1] + step;
    }
    return xxx;
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
        double area = areaEstim(x0, yy, n);
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
//
//
void run_in_cpu(double hcmin,double hcmax,int nhc,double hpmin,double hpmax,int nhp,double xmin,double xmax,int nx,double r,double tau)
{
        // sample data
        int nc;
        double *callprice = readArray<double>("callprice.txt", &nc);
        double *callstrike = readArray<double>("callstrike.txt", &nc);
        double *callopenint = readArray<double>("callopenint.txt", &nc);
        //
        int np;
        double *putprice = readArray<double>("putprice.txt", &np);
        double *putstrike = readArray<double>("putstrike.txt", &np);
        double *putopenint = readArray<double>("putopenint.txt", &np);
        //
        int nu;
        double *strike = readArray<double>("strike.txt", &nu);
        //
        // end of sample data
        //
        printf("Number of calls:   %d\n",nc);
        printf("Number of puts:    %d\n",np);
        printf("Number of strikes: %d\n",nu);
        //
		//
		//
	int mRange = nx;
	double *xRange = linspace(xmin,xmax,mRange);
	//
        NPRND *nprnd = new NPRND(callprice, callstrike, callopenint, nc,  //call data
                                 putprice, putstrike, putopenint, np,     //put  data
                                 xRange,mRange,
                                 strike, nu,
				 r, tau);                             
        //
        double time_spent = 0.0;
        clock_t begin = clock();
        //
        //
        double *hcv = linspace(hcmin, hcmax, nhc);
        double *hpv = linspace(hpmin, hpmax, nhp);
        //
	double sol[2];
        nprnd->optim_bandwidth(sol,hcv, nhc, hpv, nhp);
        //
        free(hcv);
        free(hpv);
        //
        printf("hc: %.4f\n",sol[0]);
        printf("hp: %.4f\n",sol[1]);
        printf("number of solved QP problems: %lu\n",nprnd->get_nqpprob());
        printf("number of iterations: %lu\n",nprnd->get_niterations());
        //
        clock_t end = clock();
        time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
        printf("The elapsed time is %f seconds\n", time_spent);
        //
        //
        delete nprnd;
        //
        //
        // free memory
        free(callprice);
        free(callstrike);
        free(callopenint);
        //
        free(putprice);
        free(putstrike);
        free(putopenint);
        //
        free(strike);
	free(xRange);
        //
}


