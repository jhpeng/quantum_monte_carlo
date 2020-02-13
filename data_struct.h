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
    int* sigma;
    int* sequence;
    gsl_rng* rng;
}

typedef struct placeholder{
    int leg;
    int Nsite;
    int legnth;
    int* sigma;
    int* linkv;
    int* vfirst;
    int* vlast;
}

#endif
