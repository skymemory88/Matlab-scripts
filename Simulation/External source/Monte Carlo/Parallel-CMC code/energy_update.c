#include "mex.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

void update_energy(const int* ion_types,
                   const double* ion_glande,
                   const double* ion_F,
                   const double* field,
                   const int size,
                   const double* pos,
                   const double* mom,
                   const double* staggfield,
                   const int Nions,
                   const double* inter,
                   const int new_ion,
                   const double* new_mom,
                   const double new_E_Zee,
                   const double old_E_Zee,
                   double* dE);

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    const mxArray *m_ion, *m_field, *m_lattice, *m_size, *m_inter, *m_new_ion, *m_new_mom, *m_ion_i, *m_old_E_Zee, *m_new_E_Zee, *dummy;
    double *ion_glande, *ion_F, *pos, *mom, *staggfield, *dE;
    int *ion_types;
    int i, Nions, new_ion, size;
    double old_E_Zee, new_E_Zee;
//     const int* type;
    const double *field, *inter, *new_mom;
    const mwSize* dum2;
    if(nrhs != 9)
        mexErrMsgTxt("Nine input arguments are required.");
    else if(nlhs>1)
        mexErrMsgTxt("Too many output arguments.");
    m_ion=prhs[0];
    m_field=prhs[1];
    m_lattice=prhs[2];
    m_size=prhs[3];
    m_inter=prhs[4];
    m_new_ion=prhs[5];
    m_new_mom=prhs[6];
    m_new_E_Zee=prhs[7];
    m_old_E_Zee=prhs[8];
    
    /*get number of ions*/
    dum2=mxGetDimensions(m_lattice);
    Nions=dum2[1];
    
    /*copy positions and moments*/
    pos=malloc(sizeof(double)*Nions*3);
    mom=malloc(sizeof(double)*Nions*3);
    staggfield=malloc(sizeof(double)*Nions*3);
    ion_types=malloc(sizeof(double)*Nions);
    
    
    
    for(i=0;i<Nions;++i){
        m_ion_i=mxGetCell(m_lattice,i);
        dummy=mxGetField(m_ion_i,0,"position");
        memcpy(&pos[3*i],mxGetDoubles(dummy),3*sizeof(double));
        dummy=mxGetField(m_ion_i,0,"mom");
        memcpy(&mom[3*i],mxGetDoubles(dummy),3*sizeof(double));
        dummy=mxGetField(m_ion_i,0,"staggfield");
        memcpy(&staggfield[3*i],mxGetDoubles(dummy),3*sizeof(double));   
        dummy=mxGetField(m_ion_i,0,"num");
        ion_types[i]=(int)*mxGetDoubles(dummy)-1;
       
    }
    
     ion_glande=malloc(sizeof(double)*6);
     ion_F=malloc(sizeof(double)*6);
     
      for(i=0;i<6;++i){  
        dummy=mxGetField(m_ion_i,0,"gLande");
        ion_glande[i]=*mxGetDoubles(dummy);
        dummy=mxGetField(m_ion_i,0,"F");
        ion_F[i]=*mxGetDoubles(dummy);
    }
    
    /*copy parameters for the modified ion*/
    new_ion=(int)mxGetDoubles(m_new_ion)[0]-1;
    field=mxGetDoubles(m_field);
    inter=mxGetDoubles(m_inter);
    new_mom=mxGetDoubles(m_new_mom);
    size=(int)mxGetDoubles(m_size)[0];
    old_E_Zee=(double)mxGetDoubles(m_old_E_Zee)[0];
    new_E_Zee=(double)mxGetDoubles(m_new_E_Zee)[0];
    
  
    plhs[0]=mxCreateDoubleScalar(0);
    dE=mxGetDoubles(plhs[0]);

    update_energy(ion_types, ion_glande, ion_F, field, size, pos, mom, staggfield, Nions, inter, new_ion, new_mom, new_E_Zee, old_E_Zee, dE);
    
    free(pos);
    free(mom);
    free(staggfield);
    free(ion_types);
    free(ion_glande); 
    free(ion_F);
}

void update_energy(const int* ion_types, //108
                   const double* ion_glande, //6
                   const double* ion_F, //6
                   const double* field, //3
                   const int size, //1, size=1 (L=2*size+1)
                   const double* pos, //3*108
                   const double* mom, //3*108
                   const double* staggfield, //3*108
                   const int Nions, //1, Nions=108
                   const double* inter, //201078 (9*22342)
                   const int new_ion, //1
                   const double* new_mom, //3 
                   const double new_E_Zee, //1
                   const double old_E_Zee, //1
                   double* dE) //1, output
 {
    int i, j, a, b, x, y, z, idx, maxX, maxY, maxZ;
    double  r_ij[3], staggered_field[3], sum_dipole_new[3], sum_dipole_old[3], D[9], mom1[3], mom2[3], new_mom1[3], new_mom2[3];
    double mu_b, gLande1, gLande2, H_d_old, H_d_new, F, theta_old, theta_new, H_staggfield_new, H_staggfield_old, H_corr_new, H_corr_old;
    
    // cf latmod
    maxX=2*(2*(size+0.5))+1;
    maxY=2*(2*(size+0.5))+1;
    maxZ=4*(2*(size+0.75))+1;
    
    // Bohr magneton in meV
    mu_b=5.7883817555e-2;
    
    // effective model 

    gLande1=ion_glande[ion_types[new_ion]];
    F=ion_F[ion_types[new_ion]];
    
    mom1[0]=gLande1*mom[3*new_ion];
    mom1[1]=gLande1*mom[3*new_ion+1];
    mom1[2]=gLande1*mom[3*new_ion+2];

    new_mom1[0]=gLande1*new_mom[0];
    new_mom1[1]=gLande1*new_mom[1];
    new_mom1[2]=gLande1*new_mom[2];

    staggered_field[0]=staggfield[3*new_ion];
    staggered_field[1]=staggfield[3*new_ion+1];
    staggered_field[2]=staggfield[3*new_ion+2];
    
     // angle in the XY plane
    theta_old=atan2(mom[3*new_ion+1],mom[3*new_ion]);
    theta_new=atan2(new_mom[1],new_mom[0]);
    
    sum_dipole_new[0]=0;
    sum_dipole_new[1]=0;
    sum_dipole_new[2]=0;
    sum_dipole_old[0]=0;
    sum_dipole_old[1]=0;
    sum_dipole_old[2]=0;
    
    
    for(i=0;i<Nions;++i){
        
         // effective model 
        gLande2=ion_glande[ion_types[i]];
        
        mom2[0]=gLande2*mom[3*i];
        mom2[1]=gLande2*mom[3*i+1];
        mom2[2]=gLande2*mom[3*i+2];
        
        new_mom2[0]=gLande2*new_mom[0];
        new_mom2[1]=gLande2*new_mom[1];
        new_mom2[2]=gLande2*new_mom[2];
        
        
        // sum of the dipolar interactions, in each case
        for(j=0;j<3;++j){
            r_ij[j]=pos[3*i+j]-pos[3*new_ion+j];
        }
        
        // cf latmod
        x=2*(r_ij[0]+size+0.5)+1;
        y=2*(r_ij[1]+size+0.5)+1;
        z=4*(r_ij[2]+size+0.75)+1;
        idx=((x-1)*(2*maxY+1)+(y-1))*(4*maxZ+1)+z;

        
        for(a=0;a<3;++a){
            for(b=0;b<3;++b){
                D[a*3+b]=inter[9*(idx-1)+a*3+b];
            }
        }

        if (i==new_ion){
            for(j=0;j<3;++j){
                sum_dipole_new[j]=sum_dipole_new[j]+(D[3*j]*new_mom2[0]+D[3*j+1]*new_mom2[1]+D[3*j+2]*new_mom2[2]);
            }
        }
        else{
            for(j=0;j<3;++j){
                sum_dipole_new[j]=sum_dipole_new[j]+(D[3*j]*mom2[0]+D[3*j+1]*mom2[1]+D[3*j+2]*mom2[2]);
            }
        }
          
        for(j=0;j<3;++j){
                sum_dipole_old[j]=sum_dipole_old[j]+(D[3*j]*mom2[0]+D[3*j+1]*mom2[1]+D[3*j+2]*mom2[2]);
            }

    }
    // dipolar interaction hamiltonians
    H_d_old=-(mom1[0]*sum_dipole_old[0] + mom1[1]*sum_dipole_old[1] + mom1[2]*sum_dipole_old[2]);
    H_d_new=-(new_mom1[0]*sum_dipole_new[0] + new_mom1[1]*sum_dipole_new[1] + new_mom1[2]*sum_dipole_new[2]);
    // staggered field hamiltonians
    H_staggfield_old=-mu_b*(staggered_field[0]*mom1[0]+staggered_field[1]*mom1[1]+staggered_field[2]*mom1[2]);
    H_staggfield_new=-mu_b*(staggered_field[0]*new_mom1[0]+staggered_field[1]*new_mom1[1]+staggered_field[2]*new_mom1[2]);
    // h4 anisotropy in the XY plane for Er
    H_corr_old=-F*cos(4*theta_old);
    H_corr_new=-F*cos(4*theta_new);

    *dE=(H_d_new+new_E_Zee+H_staggfield_new+H_corr_new)-(H_d_old+old_E_Zee+H_staggfield_old+H_corr_old);
    
}
