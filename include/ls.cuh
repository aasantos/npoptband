#ifndef ls_cuh
#define ls_cuh

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <cuda.h>

__device__ static const int c__0 = 0;
__device__ static const int c__1 = 1;
__device__ static const int c__2 = 2;

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MAX2(a,b) ((a) >= (b) ? (a) : (b))

#define prop 0.398942280401433


class LS{
protected:
//
//
    double *y;
    double *x;
    int nobs;
    int nvar;
    double *A;
    double *b;
    int nineq;
    double *Aeq;
    double *beq;
    double neq;
//
//
public:

__host__ __device__ LS(double *y,double *x,int nobs,int nvar,double *A,double *b,int nineq)
{
    this->y = y;
    this->x = x;
    this->nobs = nobs;
    this->nvar = nvar;
    this->A = A;
    this->b = b;
    this->nineq = nineq;
    this->Aeq = NULL;
    this->beq = NULL;
    this->neq = 0;
}

__host__ __device__ LS(double *y,double *x,int nobs,int nvar,double *A,double *b,int nineq,double *Aeq,double *beq,int neq)
{
    this->y = y;
    this->x = x;
    this->nobs = nobs;
    this->nvar = nvar;
    this->A = A;
    this->b = b;
    this->nineq = nineq;
    this->Aeq = Aeq;
    this->beq = beq;
    this->neq = neq;
}


__host__ __device__ void lsineq(double *sol)
{
    int n = nvar;
    int le = nobs;
    int me = nobs;
    int lg = nineq;
    int mg = nineq;
    double *xx = (double*)malloc(n*sizeof(double));
    double xnorm = 0.0;
    double *w = (double*)malloc(((n+1)*(mg+2) + 2*mg)*sizeof(double));
    int *jw = (int*)malloc(lg*sizeof(int));
    int mode;
    lsi_(x, y, A, b, &le, &me, &lg, &mg, &n, xx, &xnorm, w, jw, &mode);
    for(int i=0;i<n;i++) sol[i] = xx[i];
    //delete zone
    free(xx);
    free(w);
    free(jw);
}
//
//
//
__host__ __device__ void lsineqeq(double *sol)
{
    int n = nvar;
    int lc = neq;
    int mc = neq;
    int le = nobs;
    int me = nobs;
    int lg = nineq;
    int mg = nineq;
    double *xx = (double*)malloc(n*sizeof(double));
    double xnorm = 0.0;
    double *w = (double*)malloc(((n+1)*(mg+2) + 2*mg + 2*mc + me +(me*mg)*(n - mc))*sizeof(double));
    int *jw = (int*)malloc(100*lg*sizeof(int));
    int mode;
    lsei_(Aeq,beq,x,y,A,b, &lc, &mc, &le, &me, &lg, &mg, &n, xx, &xnorm, w, jw, &mode);
    for(int i=0;i<n;i++) sol[i] = xx[i];
    //delete zone
    free(xx);
    free(w);
    free(jw);
}
//
//
//
//
//
private:

    /*     COPIES A VECTOR, X, TO A VECTOR, Y, with the given increments */
    __host__ __device__   static void dcopy___(int *n_, const double *dx, int incx,
                         double *dy, int incy)
    {
        int i, n = *n_;

        if (n <= 0) return;
        if (incx == 1 && incy == 1)
            memcpy(dy, dx, sizeof(double) * ((unsigned) n));
        else if (incx == 0 && incy == 1) {
            double x = dx[0];
            for (i = 0; i < n; ++i) dy[i] = x;
        }
        else {
            for (i = 0; i < n; ++i) dy[i*incy] = dx[i*incx];
        }
    } /* dcopy___ */

    /* CONSTANT TIMES A VECTOR PLUS A VECTOR. */
    __host__ __device__ static void daxpy_sl__(int *n_, const double *da_, const double *dx,
                           int incx, double *dy, int incy)
    {
        int n = *n_, i;
        double da = *da_;

        if (n <= 0 || da == 0) return;
        for (i = 0; i < n; ++i) dy[i*incy] += da * dx[i*incx];
    }

    /* dot product dx dot dy. */
    __host__ __device__ static double ddot_sl__(int *n_, double *dx, int incx, double *dy, int incy)
    {
        int n = *n_, i;
        double sum = 0;
        if (n <= 0) return 0;
        for (i = 0; i < n; ++i) sum += dx[i*incx] * dy[i*incy];
        return (double) sum;
    }

    /* compute the L2 norm of array DX of length N, stride INCX */
    __host__ __device__ static double dnrm2___(int *n_, double *dx, int incx)
    {
        int i, n = *n_;
        double xmax = 0, scale;
        double sum = 0;
        for (i = 0; i < n; ++i) {
            double xabs = fabs(dx[incx*i]);
            if (xmax < xabs) xmax = xabs;
        }
        if (xmax == 0) return 0;
        scale = 1.0 / xmax;
        for (i = 0; i < n; ++i) {
            double xs = scale * dx[incx*i];
            sum += xs * xs;
        }
        return xmax * sqrt((double) sum);
    }


    /* apply Givens rotation */
    __host__ __device__ static void dsrot_(int n, double *dx, int incx,
                       double *dy, int incy, double *c__, double *s_)
    {
        int i;
        double c = *c__, s = *s_;

        for (i = 0; i < n; ++i) {
            double x = dx[incx*i], y = dy[incy*i];
            dx[incx*i] = c * x + s * y;
            dy[incy*i] = c * y - s * x;
        }
    }

    /* construct Givens rotation */
    __host__ __device__ static void dsrotg_(double *da, double *db, double *c, double *s)
    {
        double absa, absb, roe, scale;

        absa = fabs(*da); absb = fabs(*db);
        if (absa > absb) {
            roe = *da;
            scale = absa;
        }
        else {
            roe = *db;
            scale = absb;
        }

        if (scale != 0) {
            double r, iscale = 1 / scale;
            double tmpa = (*da) * iscale, tmpb = (*db) * iscale;
            r = (roe < 0 ? -scale : scale) * sqrt((tmpa * tmpa) + (tmpb * tmpb));
            *c = *da / r; *s = *db / r;
            *da = r;
            if (*c != 0 && fabs(*c) <= *s) *db = 1 / *c;
            else *db = *s;
        }
        else {
            *c = 1;
            *s = *da = *db = 0;
        }
    }

    __host__ __device__ static void h12_(const int *mode, int *lpivot, int *l1,
                     int *m, double *u, const int *iue, double *up,
                     double *c__, const int *ice, const int *icv, const int *ncv)
    {
        /* Initialized data */

        const double one = 1.;

        /* System generated locals */
        int u_dim1, u_offset, i__1, i__2;
        double d__1;

        /* Local variables */
        double b;
        int i__, j, i2, i3, i4;
        double cl, sm;
        int incr;
        double clinv;

        /*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 */
        /*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */
        /*     CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
        /*     HOUSEHOLDER TRANSFORMATION  Q = I + U*(U**T)/B */
        /*     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 . */
        /*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
        /*     L1,M   IF L1 <= M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
        /*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M. */
        /*            IF L1 > M THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
        /*     U(),IUE,UP */
        /*            ON ENTRY TO H1 U() STORES THE PIVOT VECTOR. */
        /*            IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS. */
        /*            ON EXIT FROM H1 U() AND UP STORE QUANTITIES DEFINING */
        /*            THE VECTOR U OF THE HOUSEHOLDER TRANSFORMATION. */
        /*            ON ENTRY TO H2 U() AND UP */
        /*            SHOULD STORE QUANTITIES PREVIOUSLY COMPUTED BY H1. */
        /*            THESE WILL NOT BE MODIFIED BY H2. */
        /*     C()    ON ENTRY TO H1 OR H2 C() STORES A MATRIX WHICH WILL BE */
        /*            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER */
        /*            TRANSFORMATION IS TO BE APPLIED. */
        /*            ON EXIT C() STORES THE SET OF TRANSFORMED VECTORS. */
        /*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
        /*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
        /*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. */
        /*            IF NCV <= 0 NO OPERATIONS WILL BE DONE ON C(). */
        /* Parameter adjustments */
        u_dim1 = *iue;
        u_offset = 1 + u_dim1;
        u -= u_offset;
        --c__;

        /* Function Body */
        if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
            goto L80;
        }
        cl = (d__1 = u[*lpivot * u_dim1 + 1], fabs(d__1));
        if (*mode == 2) {
            goto L30;
        }
        /*     ****** CONSTRUCT THE TRANSFORMATION ****** */
        i__1 = *m;
        for (j = *l1; j <= i__1; ++j) {
            sm = (d__1 = u[j * u_dim1 + 1], fabs(d__1));
            /* L10: */
            cl = MAX2(sm,cl);
        }
        if (cl <= 0.0) {
            goto L80;
        }
        clinv = one / cl;
        /* Computing 2nd power */
        d__1 = u[*lpivot * u_dim1 + 1] * clinv;
        sm = d__1 * d__1;
        i__1 = *m;
        for (j = *l1; j <= i__1; ++j) {
            /* L20: */
            /* Computing 2nd power */
            d__1 = u[j * u_dim1 + 1] * clinv;
            sm += d__1 * d__1;
        }
        cl *= sqrt(sm);
        if (u[*lpivot * u_dim1 + 1] > 0.0) {
            cl = -cl;
        }
        *up = u[*lpivot * u_dim1 + 1] - cl;
        u[*lpivot * u_dim1 + 1] = cl;
        goto L40;
        /*     ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C ****** */
    L30:
        if (cl <= 0.0) {
            goto L80;
        }
    L40:
        if (*ncv <= 0) {
            goto L80;
        }
        b = *up * u[*lpivot * u_dim1 + 1];
        if (b >= 0.0) {
            goto L80;
        }
        b = one / b;
        i2 = 1 - *icv + *ice * (*lpivot - 1);
        incr = *ice * (*l1 - *lpivot);
        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j) {
            i2 += *icv;
            i3 = i2 + incr;
            i4 = i3;
            sm = c__[i2] * *up;
            i__2 = *m;
            for (i__ = *l1; i__ <= i__2; ++i__) {
                sm += c__[i3] * u[i__ * u_dim1 + 1];
                /* L50: */
                i3 += *ice;
            }
            if (sm == 0.0) {
                goto L70;
            }
            sm *= b;
            c__[i2] += sm * *up;
            i__2 = *m;
            for (i__ = *l1; i__ <= i__2; ++i__) {
                c__[i4] += sm * u[i__ * u_dim1 + 1];
                /* L60: */
                i4 += *ice;
            }
        L70:
            ;
        }
    L80:
        return;
    } /* h12_ */



    //////////////////////////////////////////
    __host__ __device__ static void nnls_(double *a, int *mda, int *m, int *
                      n, double *b, double *x, double *rnorm, double *w,
                      double *z__, int *indx, int *mode)
    {
        /* Initialized data */

        const double one = 1.;
        const double factor = .01;

        /* System generated locals */
        int a_dim1, a_offset, i__1, i__2;
        double d__1;

        /* Local variables */
        double c__;
        int i__, j, k, l;
        double s, t;
        int ii, jj = 0, ip, iz, jz;
        double up;
        int iz1, iz2, npp1, iter;
        double wmax, alpha, asave;
        int itmax, izmax = 0, nsetp;
        double unorm;

        /*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY: */
        /*     'SOLVING LEAST SQUARES PROBLEMS'. PRENTICE-HALL.1974 */
        /*      **********   NONNEGATIVE LEAST SQUARES   ********** */
        /*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B, COMPUTE AN */
        /*     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM */
        /*                  A*X = B  SUBJECT TO  X >= 0 */
        /*     A(),MDA,M,N */
        /*            MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE ARRAY,A(). */
        /*            ON ENTRY A()  CONTAINS THE M BY N MATRIX,A. */
        /*            ON EXIT A() CONTAINS THE PRODUCT Q*A, */
        /*            WHERE Q IS AN M BY M ORTHOGONAL MATRIX GENERATED */
        /*            IMPLICITLY BY THIS SUBROUTINE. */
        /*            EITHER M>=N OR M<N IS PERMISSIBLE. */
        /*            THERE IS NO RESTRICTION ON THE RANK OF A. */
        /*     B()    ON ENTRY B() CONTAINS THE M-VECTOR, B. */
        /*            ON EXIT B() CONTAINS Q*B. */
        /*     X()    ON ENTRY X() NEED NOT BE INITIALIZED. */
        /*            ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR. */
        /*     RNORM  ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
        /*            RESIDUAL VECTOR. */
        /*     W()    AN N-ARRAY OF WORKING SPACE. */
        /*            ON EXIT W() WILL CONTAIN THE DUAL SOLUTION VECTOR. */
        /*            W WILL SATISFY W(I)=0 FOR ALL I IN SET P */
        /*            AND W(I)<=0 FOR ALL I IN SET Z */
        /*     Z()    AN M-ARRAY OF WORKING SPACE. */
        /*     INDX()AN INT WORKING ARRAY OF LENGTH AT LEAST N. */
        /*            ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
        /*            P AND Z AS FOLLOWS: */
        /*            INDX(1)    THRU INDX(NSETP) = SET P. */
        /*            INDX(IZ1)  THRU INDX (IZ2)  = SET Z. */
        /*            IZ1=NSETP + 1 = NPP1, IZ2=N. */
        /*     MODE   THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANING: */
        /*            1    THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
        /*            2    THE DIMENSIONS OF THE PROBLEM ARE WRONG, */
        /*                 EITHER M <= 0 OR N <= 0. */
        /*            3    ITERATION COUNT EXCEEDED, MORE THAN 3*N ITERATIONS. */
        /* Parameter adjustments */
        --z__;
        --b;
        --indx;
        --w;
        --x;
        a_dim1 = *mda;
        a_offset = 1 + a_dim1;
        a -= a_offset;

        /* Function Body */
        /*     revised          Dieter Kraft, March 1983 */
        *mode = 2;
        if (*m <= 0 || *n <= 0) {
            goto L290;
        }
        *mode = 1;
        iter = 0;
        itmax = *n * 3;
        /* STEP ONE (INITIALIZE) */
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* L100: */
            indx[i__] = i__;
        }
        iz1 = 1;
        iz2 = *n;
        nsetp = 0;
        npp1 = 1;
        x[1] = 0.0;
        dcopy___(n, &x[1], 0, &x[1], 1);
        /* STEP TWO (COMPUTE DUAL VARIABLES) */
        /* .....ENTRY LOOP A */
    L110:
        if (iz1 > iz2 || nsetp >= *m) {
            goto L280;
        }
        i__1 = iz2;
        for (iz = iz1; iz <= i__1; ++iz) {
            j = indx[iz];
            /* L120: */
            i__2 = *m - nsetp;
            w[j] = ddot_sl__(&i__2, &a[npp1 + j * a_dim1], 1, &b[npp1], 1)
            ;
        }
        /* STEP THREE (TEST DUAL VARIABLES) */
    L130:
        wmax = 0.0;
        i__2 = iz2;
        for (iz = iz1; iz <= i__2; ++iz) {
            j = indx[iz];
            if (w[j] <= wmax) {
                goto L140;
            }
            wmax = w[j];
            izmax = iz;
        L140:
            ;
        }
        /* .....EXIT LOOP A */
        if (wmax <= 0.0) {
            goto L280;
        }
        iz = izmax;
        j = indx[iz];
        /* STEP FOUR (TEST INDX J FOR LINEAR DEPENDENCY) */
        asave = a[npp1 + j * a_dim1];
        i__2 = npp1 + 1;
        h12_(&c__1, &npp1, &i__2, m, &a[j * a_dim1 + 1], &c__1, &up, &z__[1], &
             c__1, &c__1, &c__0);
        unorm = dnrm2___(&nsetp, &a[j * a_dim1 + 1], 1);
        t = factor * (d__1 = a[npp1 + j * a_dim1], fabs(d__1));
        d__1 = unorm + t;
        if (d__1 - unorm <= 0.0) {
            goto L150;
        }
        dcopy___(m, &b[1], 1, &z__[1], 1);
        i__2 = npp1 + 1;
        h12_(&c__2, &npp1, &i__2, m, &a[j * a_dim1 + 1], &c__1, &up, &z__[1], &
             c__1, &c__1, &c__1);
        if (z__[npp1] / a[npp1 + j * a_dim1] > 0.0) {
            goto L160;
        }
    L150:
        a[npp1 + j * a_dim1] = asave;
        w[j] = 0.0;
        goto L130;
        /* STEP FIVE (ADD COLUMN) */
    L160:
        dcopy___(m, &z__[1], 1, &b[1], 1);
        indx[iz] = indx[iz1];
        indx[iz1] = j;
        ++iz1;
        nsetp = npp1;
        ++npp1;
        i__2 = iz2;
        for (jz = iz1; jz <= i__2; ++jz) {
            jj = indx[jz];
            /* L170: */
            h12_(&c__2, &nsetp, &npp1, m, &a[j * a_dim1 + 1], &c__1, &up, &a[jj *
                                                                             a_dim1 + 1], &c__1, mda, &c__1);
        }
        k = MIN2(npp1,*mda);
        w[j] = 0.0;
        i__2 = *m - nsetp;
        dcopy___(&i__2, &w[j], 0, &a[k + j * a_dim1], 1);
        /* STEP SIX (SOLVE LEAST SQUARES SUB-PROBLEM) */
        /* .....ENTRY LOOP B */
    L180:
        for (ip = nsetp; ip >= 1; --ip) {
            if (ip == nsetp) {
                goto L190;
            }
            d__1 = -z__[ip + 1];
            daxpy_sl__(&ip, &d__1, &a[jj * a_dim1 + 1], 1, &z__[1], 1);
        L190:
            jj = indx[ip];
            /* L200: */
            z__[ip] /= a[ip + jj * a_dim1];
        }
        ++iter;
        if (iter <= itmax) {
            goto L220;
        }
    L210:
        *mode = 3;
        goto L280;
        /* STEP SEVEN TO TEN (STEP LENGTH ALGORITHM) */
    L220:
        alpha = one;
        jj = 0;
        i__2 = nsetp;
        for (ip = 1; ip <= i__2; ++ip) {
            if (z__[ip] > 0.0) {
                goto L230;
            }
            l = indx[ip];
            t = -x[l] / (z__[ip] - x[l]);
            if (alpha < t) {
                goto L230;
            }
            alpha = t;
            jj = ip;
        L230:
            ;
        }
        i__2 = nsetp;
        for (ip = 1; ip <= i__2; ++ip) {
            l = indx[ip];
            /* L240: */
            x[l] = (one - alpha) * x[l] + alpha * z__[ip];
        }
        /* .....EXIT LOOP B */
        if (jj == 0) {
            goto L110;
        }
        /* STEP ELEVEN (DELETE COLUMN) */
        i__ = indx[jj];
    L250:
        x[i__] = 0.0;
        ++jj;
        i__2 = nsetp;
        for (j = jj; j <= i__2; ++j) {
            ii = indx[j];
            indx[j - 1] = ii;
            dsrotg_(&a[j - 1 + ii * a_dim1], &a[j + ii * a_dim1], &c__, &s);
            t = a[j - 1 + ii * a_dim1];
            dsrot_(*n, &a[j - 1 + a_dim1], *mda, &a[j + a_dim1], *mda, &c__, &s);
            a[j - 1 + ii * a_dim1] = t;
            a[j + ii * a_dim1] = 0.0;
            /* L260: */
            dsrot_(1, &b[j - 1], 1, &b[j], 1, &c__, &s);
        }
        npp1 = nsetp;
        --nsetp;
        --iz1;
        indx[iz1] = i__;
        if (nsetp <= 0) {
            goto L210;
        }
        i__2 = nsetp;
        for (jj = 1; jj <= i__2; ++jj) {
            i__ = indx[jj];
            if (x[i__] <= 0.0) {
                goto L250;
            }
            /* L270: */
        }
        dcopy___(m, &b[1], 1, &z__[1], 1);
        goto L180;
        /* STEP TWELVE (SOLUTION) */
    L280:
        k = MIN2(npp1,*m);
        i__2 = *m - nsetp;
        *rnorm = dnrm2___(&i__2, &b[k], 1);
        if (npp1 > *m) {
            w[1] = 0.0;
            dcopy___(n, &w[1], 0, &w[1], 1);
        }
        /* END OF SUBROUTINE NNLS */
    L290:
        return;
    } /* nnls_ */



    ////////////////////////////////////
    __host__ __device__ static void ldp_(double *g, int *mg, int *m, int *n,
                     double *h__, double *x, double *xnorm, double *w,
                     int *indx, int *mode)
    {
        /* Initialized data */

        const double one = 1.;

        /* System generated locals */
        int g_dim1, g_offset, i__1, i__2;
        double d__1;

        /* Local variables */
        int i__, j, n1, if__, iw, iy, iz;
        double fac;
        double rnorm;
        int iwdual;

        /*                     T                                                */
        /*     MINIMIZE   1/2 X X    SUBJECT TO   G * X >= H.                   */
        /*       C.L. LAWSON, R.J. HANSON: 'SOLVING LEAST SQUARES PROBLEMS'     */
        /*       PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY, 1974.             */
        /*     PARAMETER DESCRIPTION:                                           */
        /*     G(),MG,M,N   ON ENTRY G() STORES THE M BY N MATRIX OF            */
        /*                  LINEAR INEQUALITY CONSTRAINTS. G() HAS FIRST        */
        /*                  DIMENSIONING PARAMETER MG                           */
        /*     H()          ON ENTRY H() STORES THE M VECTOR H REPRESENTING     */
        /*                  THE RIGHT SIDE OF THE INEQUALITY SYSTEM             */
        /*     REMARK: G(),H() WILL NOT BE CHANGED DURING CALCULATIONS BY LDP   */
        /*     X()          ON ENTRY X() NEED NOT BE INITIALIZED.               */
        /*                  ON EXIT X() STORES THE SOLUTION VECTOR X IF MODE=1. */
        /*     XNORM        ON EXIT XNORM STORES THE EUCLIDIAN NORM OF THE      */
        /*                  SOLUTION VECTOR IF COMPUTATION IS SUCCESSFUL        */
        /*     W()          W IS A ONE DIMENSIONAL WORKING SPACE, THE LENGTH    */
        /*                  OF WHICH SHOULD BE AT LEAST (M+2)*(N+1) + 2*M       */
        /*                  ON EXIT W() STORES THE LAGRANGE MULTIPLIERS         */
        /*                  ASSOCIATED WITH THE CONSTRAINTS                     */
        /*                  AT THE SOLUTION OF PROBLEM LDP                      */
        /*     INDX()      INDX() IS A ONE DIMENSIONAL INT WORKING SPACE        */
        /*                  OF LENGTH AT LEAST M                                */
        /*     MODE         MODE IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING   */
        /*                  MEANINGS:                                           */
        /*          MODE=1: SUCCESSFUL COMPUTATION                              */
        /*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N.LE.0)   */
        /*               3: ITERATION COUNT EXCEEDED BY NNLS */
        /*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE */
        /* Parameter adjustments */
        --indx;
        --h__;
        --x;
        g_dim1 = *mg;
        g_offset = 1 + g_dim1;
        g -= g_offset;
        --w;

        /* Function Body */
        *mode = 2;
        if (*n <= 0) {
            goto L50;
        }
        /*  STATE DUAL PROBLEM */
        *mode = 1;
        x[1] = 0.0;
        dcopy___(n, &x[1], 0, &x[1], 1);
        *xnorm = 0.0;
        if (*m == 0) {
            goto L50;
        }
        iw = 0;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                ++iw;
                /* L10: */
                w[iw] = g[j + i__ * g_dim1];
            }
            ++iw;
            /* L20: */
            w[iw] = h__[j];
        }
        if__ = iw + 1;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            ++iw;
            /* L30: */
            w[iw] = 0.0;
        }
        w[iw + 1] = one;
        n1 = *n + 1;
        iz = iw + 2;
        iy = iz + n1;
        iwdual = iy + *m;
        /*  SOLVE DUAL PROBLEM */
        nnls_(&w[1], &n1, &n1, m, &w[if__], &w[iy], &rnorm, &w[iwdual], &w[iz], &
              indx[1], mode);
        if (*mode != 1) {
            goto L50;
        }
        *mode = 4;
        if (rnorm <= 0.0) {
            goto L50;
        }
        /*  COMPUTE SOLUTION OF PRIMAL PROBLEM */
        fac = one - ddot_sl__(m, &h__[1], 1, &w[iy], 1);
        d__1 = one + fac;
        if (d__1 - one <= 0.0) {
            goto L50;
        }
        *mode = 1;
        fac = one / fac;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            /* L40: */
            x[j] = fac * ddot_sl__(m, &g[j * g_dim1 + 1], 1, &w[iy], 1);
        }
        *xnorm = dnrm2___(n, &x[1], 1);
        /*  COMPUTE LAGRANGE MULTIPLIERS FOR PRIMAL PROBLEM */
        w[1] = 0.0;
        dcopy___(m, &w[1], 0, &w[1], 1);
        daxpy_sl__(m, &fac, &w[iy], 1, &w[1], 1);
        /*  END OF SUBROUTINE LDP */
    L50:
        return;
    } /* ldp_ */


    ////////////////////////////////////////


    __host__ __device__ static void lsi_(double *e, double *f, double *g,
                     double *h__, int *le, int *me, int *lg, int *mg,
                     int *n, double *x, double *xnorm, double *w, int *
                     jw, int *mode)
    {
        /* Initialized data */

        const double epmach = 2.22e-16;
        const double one = 1.;

        /* System generated locals */
        int e_dim1, e_offset, g_dim1, g_offset, i__1, i__2, i__3;
        double d__1;

        /* Local variables */
        int i__, j;
        double t;

        /*     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF             */
        /*     INEQUALITY CONSTRAINED LINEAR LEAST SQUARES PROBLEM:             */
        /*                    MIN ||E*X-F||                                     */
        /*                     X                                                */
        /*                    S.T.  G*X >= H                                    */
        /*     THE ALGORITHM IS BASED ON QR DECOMPOSITION AS DESCRIBED IN       */
        /*     CHAPTER 23.5 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS  */
        /*     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM      */
        /*     ARE NECESSARY                                                    */
        /*     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)                       */
        /*     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )                       */
        /*     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)                       */
        /*     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )                       */
        /*     DIM(X) :   N                                                     */
        /*     DIM(W) :   (N+1)*(MG+2) + 2*MG                                   */
        /*     DIM(JW):   LG                                                    */
        /*     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS E, F, G, AND H.     */
        /*     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.         */
        /*     X     STORES THE SOLUTION VECTOR                                 */
        /*     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM      */
        /*     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST     */
        /*           MG ELEMENTS                                                */
        /*     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:     */
        /*          MODE=1: SUCCESSFUL COMPUTATION                              */
        /*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)      */
        /*               3: ITERATION COUNT EXCEEDED BY NNLS                    */
        /*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE                 */
        /*               5: MATRIX E IS NOT OF FULL RANK                        */
        /*     03.01.1980, DIETER KRAFT: CODED                                  */
        /*     20.03.1987, DIETER KRAFT: REVISED TO FORTRAN 77                  */
        /* Parameter adjustments                                                */
        --f;
        --jw;
        --h__;
        --x;
        g_dim1 = *lg;
        g_offset = 1 + g_dim1;
        g -= g_offset;
        e_dim1 = *le;
        e_offset = 1 + e_dim1;
        e -= e_offset;
        --w;

        /* Function Body */
        /*  QR-FACTORS OF E AND APPLICATION TO F */
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* Computing MIN */
            i__2 = i__ + 1;
            j = MIN2(i__2,*n);
            i__2 = i__ + 1;
            i__3 = *n - i__;
            h12_(&c__1, &i__, &i__2, me, &e[i__ * e_dim1 + 1], &c__1, &t, &e[j *
                                                                             e_dim1 + 1], &c__1, le, &i__3);
            /* L10: */
            i__2 = i__ + 1;
            h12_(&c__2, &i__, &i__2, me, &e[i__ * e_dim1 + 1], &c__1, &t, &f[1], &
                 c__1, &c__1, &c__1);
        }
        /*  TRANSFORM G AND H TO GET LEAST DISTANCE PROBLEM */
        *mode = 5;
        i__2 = *mg;
        for (i__ = 1; i__ <= i__2; ++i__) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if ((d__1 = e[j + j * e_dim1], fabs(d__1)) < epmach) {
                    goto L50;
                }
                /* L20: */
                i__3 = j - 1;
                g[i__ + j * g_dim1] = (g[i__ + j * g_dim1] - ddot_sl__(&i__3, &g[
                                                                                 i__ + g_dim1], *lg, &e[j * e_dim1 + 1], 1)) / e[j + j *
                                                                                                                                 e_dim1];
            }
            /* L30: */
            h__[i__] -= ddot_sl__(n, &g[i__ + g_dim1], *lg, &f[1], 1);
        }
        /*  SOLVE LEAST DISTANCE PROBLEM */
        ldp_(&g[g_offset], lg, mg, n, &h__[1], &x[1], xnorm, &w[1], &jw[1], mode);
        if (*mode != 1) {
            goto L50;
        }
        /*  SOLUTION OF ORIGINAL PROBLEM */
        daxpy_sl__(n, &one, &f[1], 1, &x[1], 1);
        for (i__ = *n; i__ >= 1; --i__) {
            /* Computing MIN */
            i__2 = i__ + 1;
            j = MIN2(i__2,*n);
            /* L40: */
            i__2 = *n - i__;
            x[i__] = (x[i__] - ddot_sl__(&i__2, &e[i__ + j * e_dim1], *le, &x[j], 1))
            / e[i__ + i__ * e_dim1];
        }
        /* Computing MIN */
        i__2 = *n + 1;
        j = MIN2(i__2,*me);
        i__2 = *me - *n;
        t = dnrm2___(&i__2, &f[j], 1);
        *xnorm = sqrt(*xnorm * *xnorm + t * t);
        /*  END OF SUBROUTINE LSI */
    L50:
        return;
    } /* lsi_ */



    ////////////////////
    __host__ __device__ static void hfti_(double *a, int *mda, int *m, int *
                      n, double *b, int *mdb, const int *nb, double *tau, int
                      *krank, double *rnorm, double *h__, double *g, int *
                      ip)
    {
        /* Initialized data */

        const double factor = .001;

        /* System generated locals */
        int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
        double d__1;

        /* Local variables */
        int i__, j, k, l;
        int jb, kp1;
        double tmp, hmax = 0.0;
        int lmax, ldiag;

        /*     RANK-DEFICIENT LEAST SQUARES ALGORITHM AS DESCRIBED IN: */
        /*     C.L.LAWSON AND R.J.HANSON, JET PROPULSION LABORATORY, 1973 JUN 12 */
        /*     TO APPEAR IN 'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974 */
        /*     A(*,*),MDA,M,N   THE ARRAY A INITIALLY CONTAINS THE M x N MATRIX A */
        /*                      OF THE LEAST SQUARES PROBLEM AX = B. */
        /*                      THE FIRST DIMENSIONING PARAMETER MDA MUST SATISFY */
        /*                      MDA >= M. EITHER M >= N OR M < N IS PERMITTED. */
        /*                      THERE IS NO RESTRICTION ON THE RANK OF A. */
        /*                      THE MATRIX A WILL BE MODIFIED BY THE SUBROUTINE. */
        /*     B(*,*),MDB,NB    IF NB = 0 THE SUBROUTINE WILL MAKE NO REFERENCE */
        /*                      TO THE ARRAY B. IF NB > 0 THE ARRAY B() MUST */
        /*                      INITIALLY CONTAIN THE M x NB MATRIX B  OF THE */
        /*                      THE LEAST SQUARES PROBLEM AX = B AND ON RETURN */
        /*                      THE ARRAY B() WILL CONTAIN THE N x NB SOLUTION X. */
        /*                      IF NB>1 THE ARRAY B() MUST BE DOUBLE SUBSCRIPTED */
        /*                      WITH FIRST DIMENSIONING PARAMETER MDB>=MAX(M,N), */
        /*                      IF NB=1 THE ARRAY B() MAY BE EITHER SINGLE OR */
        /*                      DOUBLE SUBSCRIPTED. */
        /*     TAU              ABSOLUTE TOLERANCE PARAMETER FOR PSEUDORANK */
        /*                      DETERMINATION, PROVIDED BY THE USER. */
        /*     KRANK            PSEUDORANK OF A, SET BY THE SUBROUTINE. */
        /*     RNORM            ON EXIT, RNORM(J) WILL CONTAIN THE EUCLIDIAN */
        /*                      NORM OF THE RESIDUAL VECTOR FOR THE PROBLEM */
        /*                      DEFINED BY THE J-TH COLUMN VECTOR OF THE ARRAY B. */
        /*     H(), G()         ARRAYS OF WORKING SPACE OF LENGTH >= N. */
        /*     IP()             INT ARRAY OF WORKING SPACE OF LENGTH >= N */
        /*                      RECORDING PERMUTATION INDICES OF COLUMN VECTORS */
        /* Parameter adjustments */
        --ip;
        --g;
        --h__;
        a_dim1 = *mda;
        a_offset = 1 + a_dim1;
        a -= a_offset;
        --rnorm;
        b_dim1 = *mdb;
        b_offset = 1 + b_dim1;
        b -= b_offset;

        /* Function Body */
        k = 0;
        ldiag = MIN2(*m,*n);
        if (ldiag <= 0) {
            goto L270;
        }
        /*   COMPUTE LMAX */
        i__1 = ldiag;
        for (j = 1; j <= i__1; ++j) {
            if (j == 1) {
                goto L20;
            }
            lmax = j;
            i__2 = *n;
            for (l = j; l <= i__2; ++l) {
                /* Computing 2nd power */
                d__1 = a[j - 1 + l * a_dim1];
                h__[l] -= d__1 * d__1;
                /* L10: */
                if (h__[l] > h__[lmax]) {
                    lmax = l;
                }
            }
            d__1 = hmax + factor * h__[lmax];
            if (d__1 - hmax > 0.0) {
                goto L50;
            }
        L20:
            lmax = j;
            i__2 = *n;
            for (l = j; l <= i__2; ++l) {
                h__[l] = 0.0;
                i__3 = *m;
                for (i__ = j; i__ <= i__3; ++i__) {
                    /* L30: */
                    /* Computing 2nd power */
                    d__1 = a[i__ + l * a_dim1];
                    h__[l] += d__1 * d__1;
                }
                /* L40: */
                if (h__[l] > h__[lmax]) {
                    lmax = l;
                }
            }
            hmax = h__[lmax];
            /*   COLUMN INTERCHANGES IF NEEDED */
        L50:
            ip[j] = lmax;
            if (ip[j] == j) {
                goto L70;
            }
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                tmp = a[i__ + j * a_dim1];
                a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
                /* L60: */
                a[i__ + lmax * a_dim1] = tmp;
            }
            h__[lmax] = h__[j];
            /*   J-TH TRANSFORMATION AND APPLICATION TO A AND B */
        L70:
            /* Computing MIN */
            i__2 = j + 1;
            i__ = MIN2(i__2,*n);
            i__2 = j + 1;
            i__3 = *n - j;
            h12_(&c__1, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &a[i__ *
                                                                             a_dim1 + 1], &c__1, mda, &i__3);
            /* L80: */
            i__2 = j + 1;
            h12_(&c__2, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &b[
                                                                             b_offset], &c__1, mdb, nb);
        }
        /*   DETERMINE PSEUDORANK */
        i__2 = ldiag;
        for (j = 1; j <= i__2; ++j) {
            /* L90: */
            if ((d__1 = a[j + j * a_dim1], fabs(d__1)) <= *tau) {
                goto L100;
            }
        }
        k = ldiag;
        goto L110;
    L100:
        k = j - 1;
    L110:
        kp1 = k + 1;
        /*   NORM OF RESIDUALS */
        i__2 = *nb;
        for (jb = 1; jb <= i__2; ++jb) {
            /* L130: */
            i__1 = *m - k;
            rnorm[jb] = dnrm2___(&i__1, &b[kp1 + jb * b_dim1], 1);
        }
        if (k > 0) {
            goto L160;
        }
        i__1 = *nb;
        for (jb = 1; jb <= i__1; ++jb) {
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L150: */
                b[i__ + jb * b_dim1] = 0.0;
            }
        }
        goto L270;
    L160:
        if (k == *n) {
            goto L180;
        }
        /*   HOUSEHOLDER DECOMPOSITION OF FIRST K ROWS */
        for (i__ = k; i__ >= 1; --i__) {
            /* L170: */
            i__2 = i__ - 1;
            h12_(&c__1, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &a[
                                                                          a_offset], mda, &c__1, &i__2);
        }
    L180:
        i__2 = *nb;
        for (jb = 1; jb <= i__2; ++jb) {
            /*   SOLVE K*K TRIANGULAR SYSTEM */
            for (i__ = k; i__ >= 1; --i__) {
                /* Computing MIN */
                i__1 = i__ + 1;
                j = MIN2(i__1,*n);
                /* L210: */
                i__1 = k - i__;
                b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1] - ddot_sl__(&i__1, &
                                                                         a[i__ + j * a_dim1], *mda, &b[j + jb * b_dim1], 1)) /
                a[i__ + i__ * a_dim1];
            }
            /*   COMPLETE SOLUTION VECTOR */
            if (k == *n) {
                goto L240;
            }
            i__1 = *n;
            for (j = kp1; j <= i__1; ++j) {
                /* L220: */
                b[j + jb * b_dim1] = 0.0;
            }
            i__1 = k;
            for (i__ = 1; i__ <= i__1; ++i__) {
                /* L230: */
                h12_(&c__2, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &b[jb *
                                                                              b_dim1 + 1], &c__1, mdb, &c__1);
            }
            /*   REORDER SOLUTION ACCORDING TO PREVIOUS COLUMN INTERCHANGES */
        L240:
            for (j = ldiag; j >= 1; --j) {
                if (ip[j] == j) {
                    goto L250;
                }
                l = ip[j];
                tmp = b[l + jb * b_dim1];
                b[l + jb * b_dim1] = b[j + jb * b_dim1];
                b[j + jb * b_dim1] = tmp;
            L250:
                ;
            }
        }
    L270:
        *krank = k;
    } /* hfti_ */



    /////////////////////////////////////////
    __host__ __device__ static void lsei_(double *c__, double *d__, double *e,
                      double *f, double *g, double *h__, int *lc, int *
                      mc, int *le, int *me, int *lg, int *mg, int *n,
                      double *x, double *xnrm, double *w, int *jw, int *
                      mode)
    {
        /* Initialized data */

        const double epmach = 2.22e-16;

        /* System generated locals */
        int c_dim1, c_offset, e_dim1, e_offset, g_dim1, g_offset, i__1, i__2,
        i__3;
        double d__1;

        /* Local variables */
        int i__, j, k, l;
        double t;
        int ie, if__, ig, iw, mc1;
        int krank;
        /*     FOR MODE=1, THE SUBROUTINE RETURNS THE SOLUTION X OF               */
        /*     EQUALITY & INEQUALITY CONSTRAINED LEAST SQUARES PROBLEM LSEI :     */
        /*                MIN ||E*X - F||                                         */
        /*                 X                                                      */
        /*                S.T.  C*X  = D,                                         */
        /*                      G*X >= H.                                         */
        /*     USING QR DECOMPOSITION & ORTHOGONAL BASIS OF NULLSPACE OF C        */
        /*     CHAPTER 23.6 OF LAWSON & HANSON: SOLVING LEAST SQUARES PROBLEMS.   */
        /*     THE FOLLOWING DIMENSIONS OF THE ARRAYS DEFINING THE PROBLEM        */
        /*     ARE NECESSARY                                                      */
        /*     DIM(E) :   FORMAL (LE,N),    ACTUAL (ME,N)                         */
        /*     DIM(F) :   FORMAL (LE  ),    ACTUAL (ME  )                         */
        /*     DIM(C) :   FORMAL (LC,N),    ACTUAL (MC,N)                         */
        /*     DIM(D) :   FORMAL (LC  ),    ACTUAL (MC  )                         */
        /*     DIM(G) :   FORMAL (LG,N),    ACTUAL (MG,N)                         */
        /*     DIM(H) :   FORMAL (LG  ),    ACTUAL (MG  )                         */
        /*     DIM(X) :   FORMAL (N   ),    ACTUAL (N   )                         */
        /*     DIM(W) :   2*MC+ME+(ME+MG)*(N-MC)  for LSEI                        */
        /*              +(N-MC+1)*(MG+2)+2*MG     for LSI                         */
        /*     DIM(JW):   MAX(MG,L)                                               */
        /*     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS C, D, E, F, G, AND H. */
        /*     ON RETURN, ALL ARRAYS WILL BE CHANGED BY THE SUBROUTINE.           */
        /*     X     STORES THE SOLUTION VECTOR                                   */
        /*     XNORM STORES THE RESIDUUM OF THE SOLUTION IN EUCLIDIAN NORM        */
        /*     W     STORES THE VECTOR OF LAGRANGE MULTIPLIERS IN ITS FIRST       */
        /*           MC+MG ELEMENTS                                               */
        /*     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:       */
        /*          MODE=1: SUCCESSFUL COMPUTATION                                */
        /*               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)        */
        /*               3: ITERATION COUNT EXCEEDED BY NNLS                      */
        /*               4: INEQUALITY CONSTRAINTS INCOMPATIBLE                   */
        /*               5: MATRIX E IS NOT OF FULL RANK                          */
        /*               6: MATRIX C IS NOT OF FULL RANK                          */
        /*               7: RANK DEFECT IN HFTI                                   */
        /*     18.5.1981, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN                    */
        /*     20.3.1987, DIETER KRAFT, DFVLR OBERPFAFFENHOFEN                    */
        /* Parameter adjustments                                                  */
        --d__;
        --f;
        --h__;
        --x;
        g_dim1 = *lg;
        g_offset = 1 + g_dim1;
        g -= g_offset;
        e_dim1 = *le;
        e_offset = 1 + e_dim1;
        e -= e_offset;
        c_dim1 = *lc;
        c_offset = 1 + c_dim1;
        c__ -= c_offset;
        --w;
        --jw;

        /* Function Body */
        *mode = 2;
        if (*mc > *n) {
            goto L75;
        }
        l = *n - *mc;
        mc1 = *mc + 1;
        iw = (l + 1) * (*mg + 2) + (*mg << 1) + *mc;
        ie = iw + *mc + 1;
        if__ = ie + *me * l;
        ig = if__ + *me;
        /*  TRIANGULARIZE C AND APPLY FACTORS TO E AND G */
        i__1 = *mc;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* Computing MIN */
            i__2 = i__ + 1;
            j = MIN2(i__2,*lc);
            i__2 = i__ + 1;
            i__3 = *mc - i__;
            h12_(&c__1, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &
                 c__[j + c_dim1], lc, &c__1, &i__3);
            i__2 = i__ + 1;
            h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &e[
                                                                                 e_offset], le, &c__1, me);
            /* L10: */
            i__2 = i__ + 1;
            h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &g[
                                                                                 g_offset], lg, &c__1, mg);
        }
        /*  SOLVE C*X=D AND MODIFY F */
        *mode = 6;
        i__2 = *mc;
        for (i__ = 1; i__ <= i__2; ++i__) {
            if ((d__1 = c__[i__ + i__ * c_dim1], fabs(d__1)) < epmach) {
                goto L75;
            }
            i__1 = i__ - 1;
            x[i__] = (d__[i__] - ddot_sl__(&i__1, &c__[i__ + c_dim1], *lc, &x[1], 1))
            / c__[i__ + i__ * c_dim1];
            /* L15: */
        }
        *mode = 1;
        w[mc1] = 0.0;
        i__2 = *mg; /* BUGFIX for *mc == *n: changed from *mg - *mc, SGJ 2010 */
        dcopy___(&i__2, &w[mc1], 0, &w[mc1], 1);
        if (*mc == *n) {
            goto L50;
        }
        i__2 = *me;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L20: */
            w[if__ - 1 + i__] = f[i__] - ddot_sl__(mc, &e[i__ + e_dim1], *le, &x[1], 1);
        }
        /*  STORE TRANSFORMED E & G */
        i__2 = *me;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L25: */
            dcopy___(&l, &e[i__ + mc1 * e_dim1], *le, &w[ie - 1 + i__], *me);
        }
        i__2 = *mg;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L30: */
            dcopy___(&l, &g[i__ + mc1 * g_dim1], *lg, &w[ig - 1 + i__], *mg);
        }
        if (*mg > 0) {
            goto L40;
        }
        /*  SOLVE LS WITHOUT INEQUALITY CONSTRAINTS */
        *mode = 7;
        k = MAX2(*le,*n);
        t = sqrt(epmach);
        hfti_(&w[ie], me, me, &l, &w[if__], &k, &c__1, &t, &krank, xnrm, &w[1], &
              w[l + 1], &jw[1]);
        dcopy___(&l, &w[if__], 1, &x[mc1], 1);
        if (krank != l) {
            goto L75;
        }
        *mode = 1;
        goto L50;
        /*  MODIFY H AND SOLVE INEQUALITY CONSTRAINED LS PROBLEM */
    L40:
        i__2 = *mg;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L45: */
            h__[i__] -= ddot_sl__(mc, &g[i__ + g_dim1], *lg, &x[1], 1);
        }
        lsi_(&w[ie], &w[if__], &w[ig], &h__[1], me, me, mg, mg, &l, &x[mc1], xnrm,
             &w[mc1], &jw[1], mode);
        if (*mc == 0) {
            goto L75;
        }
        t = dnrm2___(mc, &x[1], 1);
        *xnrm = sqrt(*xnrm * *xnrm + t * t);
        if (*mode != 1) {
            goto L75;
        }
        /*  SOLUTION OF ORIGINAL PROBLEM AND LAGRANGE MULTIPLIERS */
    L50:
        i__2 = *me;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L55: */
            f[i__] = ddot_sl__(n, &e[i__ + e_dim1], *le, &x[1], 1) - f[i__];
        }
        i__2 = *mc;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L60: */
            d__[i__] = ddot_sl__(me, &e[i__ * e_dim1 + 1], 1, &f[1], 1) -
            ddot_sl__(mg, &g[i__ * g_dim1 + 1], 1, &w[mc1], 1);
        }
        for (i__ = *mc; i__ >= 1; --i__) {
            /* L65: */
            i__2 = i__ + 1;
            h12_(&c__2, &i__, &i__2, n, &c__[i__ + c_dim1], lc, &w[iw + i__], &x[
                                                                                 1], &c__1, &c__1, &c__1);
        }
        for (i__ = *mc; i__ >= 1; --i__) {
            /* Computing MIN */
            i__2 = i__ + 1;
            j = MIN2(i__2,*lc);
            i__2 = *mc - i__;
            w[i__] = (d__[i__] - ddot_sl__(&i__2, &c__[j + i__ * c_dim1], 1, &
                                           w[j], 1)) / c__[i__ + i__ * c_dim1];
            /* L70: */
        }
        /*  END OF SUBROUTINE LSEI */
    L75:
        return;
    } /* lsei_ */



};  // end of the LS object


#endif