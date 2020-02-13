#include <gsl/gsl_rng.h>

#ifndef data_struct.h
#define data_struct.h

typedef struct lattice_profile{
    int leg;
    int Nsite;
    int Nb;
    double beta;
    int* bond2index;
    double* bondst
} lattice_profile;

typedef struct system_state{
    int leg;
    int Nsite;
    int length;
    int noo;
    int* sigma;
    int* sequence;
    gsl_rng* rng;
} system_state;

typedef struct placeholder{
    int leg;
    int Nsite;
    int legnth;
    int* sigma;
    int* linkv;
    int* vfirst;
    int* vlast;
} placeholder;

#endif
