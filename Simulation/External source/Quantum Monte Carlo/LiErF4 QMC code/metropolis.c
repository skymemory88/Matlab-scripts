#include "mex.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


void metropolis(const double dE,
                const double T,
                double* acc);

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    const mxArray *m_dE, *m_T;
    double dE, T;
    double *acc;
    
    if(nrhs != 2)
        mexErrMsgTxt("Two input arguments are required.");
    else if(nlhs>1)
        mexErrMsgTxt("Too many output arguments.");
    
    m_dE=prhs[0];
    m_T=prhs[1];
    
    dE=mxGetPr(m_dE)[0];
    T=mxGetPr(m_T)[0];
    
    
    plhs[0]=mxCreateDoubleScalar(0);
    acc=mxGetPr(plhs[0]);


    metropolis(dE, T, acc);
}



void metropolis(const double dE,
                const double T,
                double* acc)
{
    double alpha;
    
    if(T==0){
        if(dE<0){
         *acc=1;
        }
        else{
         *acc=0;
        }
    }

    else{
        if (exp(-11.6*(dE)/T)>=1){ // 1/kB=11.6 K/meV
            *acc=1;
        }
        else{
            alpha=(double) rand() / (RAND_MAX) ;
            if (exp(-11.6*(dE)/T)>alpha){
                *acc=1;
            }
            else{
                *acc=0;
            }
        }
    }

}
