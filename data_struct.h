#include <gsl/gsl_rng.h>

#ifndef data_struct_h
#define data_struct_h

typedef struct lattice_profile{
    int leg;
    int Nsite;
    int Nb;
    double beta;
    int* bond2index;
    double* bondst;
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
    int length;
    int* sigma;
    int* linkv;
    int* vfirst;
    int* vlast;
} placeholder;

lattice_profile* create_lattice_profile(
            int leg, 
            int Nsite, 
            int Nb);

void destroy_lattice_profile(
            lattice_profile* lap);

void lattice_profile_set_beta(
            lattice_profile* lap, 
            double beta);

system_state* create_system_state(
            int leg, 
            int Nsite, 
            int length);

void destroy_system_state(
            system_state* state);

void system_state_set_rng(
            system_state* state, 
            int seed);

placeholder* create_placeholder(
            int leg, 
            int Nsite, 
            int length);

void destroy_placeholder(
            placeholder* ws);

void sequence_doubling(
            lattice_profile* lap,
            system_state* state, 
            placeholder* ws);

#endif
