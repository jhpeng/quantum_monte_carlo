#include <stdio.h>
#include <stdlib.h>
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
    Data = (double*)malloc(sizeof(double)*n*2);
    M=n;

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

int
main (void)
{
  int i,seed=1289731;
  const int n = 1024*1024*4;
  gap_estimator_setup_workspace(n);
  double* data = (double*)malloc(sizeof(double)*2*n);
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,seed);

  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;

  for (i = 0; i < n; i++)
    {
      REAL(data,i) = gsl_rng_uniform_pos(rng);
      IMAG(data,i) = gsl_rng_uniform_pos(rng);
    }


  for (i = 0; i < 0; i++)
    {
      printf ("%d: %e %e\n", i, REAL(data,i),
                                IMAG(data,i));
    }
  printf ("\n");

  wavetable = gsl_fft_complex_wavetable_alloc (n);
  workspace = gsl_fft_complex_workspace_alloc (n);

  for (i = 0; i < (int) wavetable->nf; i++)
    {
       printf ("# factor %d: %zu\n", i,
               wavetable->factor[i]);
    }

  gsl_fft_complex_forward (data, 1, n,
                           wavetable, workspace);

  for (i = 0; i < 0; i++)
    {
      printf ("%d: %e %e\n", i, REAL(data,i),
                                IMAG(data,i));
    }

  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
  free(data);
  return 0;
}
