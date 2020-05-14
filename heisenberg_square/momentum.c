#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int Init=0;
static int Nk;
static int Ns;
static double *Struct_factor;
static int N_estimate;
static double *Data;

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

void propagate_state(int* sigmap, const int* bond2index, int sp);

void momentum_setup_workspace(int nx, int ny, int nz, int nk, double* wave_vector, int n_sample){
    if(Init){
        free(Struct_factor);
        free(Data);
    }

    Struct_factor = (double*)malloc(sizeof(double)*nx*ny*nz*nk*2);
    Data = (double*)malloc(sizeof(double)*nk*2);

    int iq,ix,iy,iz;
    double qx,qy,qz;
    int n0 = nx*ny*nz;
    int n1 = nx*ny;
    for(iq=0;iq<nk;++iq){
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
    for(iq=0;iq<2*nk;++iq) Data[iq]=0;

    Init=1;
    Nk = nk;
    Ns = n0;
    N_estimate=0;
}

void momentum_free_memory(){
    free(Struct_factor);
    free(Data);
}

void momentum_collect_data(int* sequence, int length, int* sigma0, int* sigmap, int nsite, int* bond2index, double beta, int i_sample){
    int i,j,i_bond,p,sp,type,iq;
    int noo=0;
    double ms[2*Nk],msx[2*Nk],ms2[Nk];

    for(i=0;i<nsite;++i) sigmap[i]=sigma0[i];
    for(iq=0;iq<Nk;++iq){
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
            j = bond2index[i_bond*4+1];

            for(iq=0;iq<Nk;++iq){
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
            ++noo;
        }
    }

    for(iq=0;iq<Nk;++iq){
        if(noo!=0){
            REAL(msx,iq) = beta*(REAL(msx,iq)*REAL(msx,iq)+
                IMAG(msx,iq)*IMAG(msx,iq)+ms2[iq])/noo/(noo+1)/nsite/nsite*0.25;
            ms2[iq] = ms2[iq]/noo/nsite/nsite*0.25;
        }
        else{
            ms2[iq] = (REAL(ms,iq)*REAL(ms,iq)+IMAG(ms,iq)*IMAG(ms,iq))/nsite/nsite*0.25;
            REAL(msx,iq) = beta*ms2[iq];
        }

        Data[2*iq+0] += ms2[iq];
        Data[2*iq+1] += REAL(msx,iq);
    }

    N_estimate++;
}

void momentum_calc_mean_fileout(char* prefix){
    int i;

    for(i=0;i<Nk;++i){
        Data[2*i+0] = Data[2*i+0]/N_estimate;
        Data[2*i+1] = Data[2*i+1]/N_estimate;
    }

    char filename[128];
    sprintf(filename,"%s.mom",prefix);
    FILE* file = fopen(filename,"a");
    for(i=0;i<2*Nk;++i) fprintf(file,"%.16e ",Data[i]);
    fprintf(file,"\n");

    fclose(file);
    for(i=0;i<2*Nk;++i) Data[i]=0;

    N_estimate=0;
}
