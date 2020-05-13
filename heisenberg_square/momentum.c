#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int Init=0;
static int Nq;
static int Ns;
static double *Struct_factor;
static int N_estimate;
static double *Data;

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

void propagate_state(int* sigmap, const int* bond2index, int sp);

void momentum_setup_workspace(int nx, int ny, int nz, int nq, double* wave_vector, int n_sample){
    if(Init){
        free(Struct_factor);
        free(Data);
    }

    Struct_factor = (double*)malloc(sizeof(double)*nx*ny*nz*nq*2);
    Data = (double*)malloc(sizeof(double)*n_sample*nq*2);

    int iq,ix,iy,iz;
    double qx,qy,qz;
    double n0 = nx*ny*nz;
    double n1 = nx*ny;
    for(iq=0;iq<nq;++iq){
        qx = wave_vector[3*iq+0];
        qy = wave_vector[3*iq+1];
        qz = wave_vector[3*iq+2];

        for(iz=0;iz<nz;++iz){
        for(iy=0;iy<ny;++iy){
        for(ix=0;ix<nx;++ix){
            REAL(Struct_factor,iq*n0+iz*n1+iy*nx+ix) = cos(ix*qx+iy*qy+iz*qz);
            IMAG(Struct_factor,iq*n0+iz*n1+iy*nx+ix) = sin(ix*qx+iy*qy+iz*qz);
        }
        }
        }
    }

    Init=1;
    Nq = nq;
    Ns = n0;
    N_estimate=n_sample;
}

void momentum_free_memory(){
    free(Struct_factor);
    free(Data);
}

void momentum_collect_data(int* sequence, int length, int* sigma0, int* sigmap, int nsite, int* bond2index){
    int i,j,i_bond,p,sp,type,iq;
    double ms[2*Nq],msx[2*Nq],ms2[Nq];

    for(i=0;i<nsite;++i) sigmap[i]=sigma0[i];
    for(iq=0;iq<Nq;++iq){
        REAL(ms,iq) = 0;
        IMAG(ms,iq) = 0;
        REAL(msx,iq) = 0;
        IMAG(msx,iq) = 0;
        ms2[iq] = 0;
        for(i=0;i<nsite;++i){
            REAL(ms,iq)+=sigma0[i]*REAL(Struct_factor,iq*nsite+i);
            IMAG(ms,iq)+=sigma0[i]*IMAG(Struct_factor,iq*nsite+i);
        }
    }

    for(p=0;p<length;++p){
        sp = sequence[p];

        if(sp!=-1){
            type = sp%6;
            i_bond = sp/6;
            i = bond2index[i_bond*4+0];
            j = bond2index[i_bond*4+q];

            for(iq=0;iq<Nq;++iq){
                if(type==1){
                    REAL(ms,iq) += -2*sigmap[i]*REAL(Struct_factor,iq*nsite+i);
                    REAL(ms,iq) += -2*sigmap[j]*REAL(Struct_factor,iq*nsite+j);
                    IMAG(ms,iq) += -2*sigmap[i]*IMAG(Struct_factor,iq*nsite+i);
                    IMAG(ms,iq) += -2*sigmap[j]*IMAG(Struct_factor,iq*nsite+j);
                }

                REAL(msx,iq) += REAL(ms,iq);
                IMAG(msx,iq) += IMAG(ms,iq);
                ms2[iq] += REAL(ms,iq)*REAL(ms,iq)+IMAG(ms,iq)*IMAG(ms,iq);
            }

            propagate_state(sigmap,bond2index,sp);
        }
    }
}

