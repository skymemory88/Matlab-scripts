#include "mex.h"
#include <math.h>

void dipole_sum(const double *rij, const double *L, const double *N, double* D);

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    const mxArray *rij, *L, *N;
    double *drij, *dL, *dN, *D;
    /* check number of input and output arguments*/
    if(nrhs != 3)
        mexErrMsgTxt("Three input arguments are required.");
    else if(nlhs>1)
        mexErrMsgTxt("Too many output arguments.");
    
    rij=prhs[0];
    L=prhs[1];
    N=prhs[2];
    
    /* check dimensions of input arguments */
    if(mxGetNumberOfDimensions(rij)>2 || mxGetM(rij)!=3 || mxGetN(rij)!=1)
        mexErrMsgTxt("Argument 1 must be a 3x1 vector.");
    else if(mxGetNumberOfDimensions(L)>2 || mxGetM(L)!=3 || mxGetN(L)!=1)
        mexErrMsgTxt("Argument 2 must be a 3x1 vector.");
    else if(mxGetNumberOfDimensions(N)>2 || mxGetM(N)!=3 || mxGetN(N)!=1)
        mexErrMsgTxt("Argument 3 must be a 3x1 vector.");
    
    plhs[0]=mxCreateNumericMatrix(3,3,mxDOUBLE_CLASS,mxREAL);
    drij=mxGetDoubles(rij);
    dL=mxGetDoubles(L);
    dN=mxGetDoubles(N);
    D=mxGetDoubles(plhs[0]);
    dipole_sum(drij,dL,dN,D);
}

void dipole_sum(const double *rij, const double *L, const double *N, double *D)
{
    double R[3];
    double Rn;
    double V=L[0]*L[1]*L[2];
    int h,k,l,a,b;
    for(h=-N[0];h<=N[0];++h){
        for(k=-N[1];k<=N[1];++k){
            for(l=-N[2];l<=N[2];++l){
                // printf("h = %d, k = %d, l = %d \n",h,k,l);
                R[0]=rij[0]+h*L[0];
                R[1]=rij[1]+k*L[1];
                R[2]=rij[2]+l*L[2];
                Rn=sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
                // arg1: 1e-4 appears to be the lower cut-off distance for dipole calculations --Yikai
                // arg2: Pseudo-open boundary condition (upper cut-off distance) --Yikai
                if(Rn>1e-4 && Rn<N[0]*L[0]){
                    for(a=0;a<3;++a){
                        for(b=0;b<3;++b){
                            D[a*3+b]+=3*R[a]*R[b]/pow(Rn, 5)-(a==b)/pow(Rn, 3);
                        }
                    }
                }
            }
        }
    }
    /*Lorentz factor on the diagonal, finally added later*/
    /*D[0]+=4*M_PI/3/V;
    D[4]+=4*M_PI/3/V;
    D[8]+=4*M_PI/3/V;*/
}