#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_rng.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

static int Initialize=0;
static gsl_fft_complex_wavetable* Wavetable;
static gsl_fft_complex_workspace* Workspace;
static double* Powerspectrum;
static double* Data;
static int M;
static int N_estimate;

void propagate_state(int* sigmap, const int* bond2index, int sp);

void gap_estimator_setup_workspace(int length){
    int npow = (int)log2((double)length);
    int n = (int)pow(2,npow-2);

    if(Initialize){
        gsl_fft_complex_wavetable_free(Wavetable);
        gsl_fft_complex_workspace_free(Workspace);
        free(Powerspectrum);
        free(Data);
    }

    Wavetable = gsl_fft_complex_wavetable_alloc(n);
    Workspace = gsl_fft_complex_workspace_alloc(n);

    Powerspectrum = (double*)malloc(sizeof(double)*n);
    for(int i=0;i<n;++i) Powerspectrum[i]=0;
    Data = (double*)malloc(sizeof(double)*n*2);
    M=n;
    N_estimate=0;

    Initialize=1;
    
}

void gap_estimator_free_memory(){
    gsl_fft_complex_wavetable_free(Wavetable);
    gsl_fft_complex_workspace_free(Workspace);
    free(Powerspectrum);
    free(Data);
}

void gap_estimator_collect_data(int* sequence, int length, int* sigma0, int* sigmap, int nsite, int* bond2index, double* structfactor){
    int i,p,sp,type;
    int k=0;
    double dist = (double)length/M;
    double ms=0;

    for(i=0;i<nsite;++i) sigmap[i]=sigma0[i];
    for(i=0;i<nsite;++i) ms+=sigma0[i]*structfactor[i];

    for(p=0;p<length;++p){
        sp = sequence[p];
        type = sp%6;

        if(type==1) ms += -4*sigmap[i]*structfactor[i];

        if(p>k*dist && k<M){
            REAL(Data,k) = ms;
            IMAG(Data,k) = 0;
            k++;
        }

        propagate_state(sigmap,bond2index,sp);
    }
}

void gap_estimator_calc_power_spectrum(){
    int i;
    double v;

    gsl_fft_complex_forward(Data,1,M,Wavetable,Workspace);

    for(i=0;i<M;++i){
        v = REAL(Data,i)*REAL(Data,i)+IMAG(Data,i)*IMAG(Data,i);
        Powerspectrum[i] += v;
    }

    N_estimate++;
}

void gap_estimator_ave_power_spectrum_fileout(char* prefix){
    int i;

    for(i=0;i<M;++i) Powerspectrum[i] = Powerspectrum[i]/N_estimate/sqrt(M);

    char filename[128];
    sprintf(filename,"%s.pow",prefix);
    FILE *file = fopen(filename,"a");

    for(int i=0;i<M;++i) fprintf(file,"%.16e ",Powerspectrum[i]);
    fprintf(file,"\n");

    for(i=0;i<M;++i) Powerspectrum[i] = 0;
    N_estimate=0;
}
