#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "data_struct.h"

lattice_profile* create_lattice_profile(int leg, int Nsite, int Nb){
    if(leg%2==1){
        printf("data_struct.c create_lattice_profile : leg should be even number\n");
        exit(-1);
    }

    int nl=leg/2;

    lattice_profile* lap = (lattice_profile*)malloc(sizeof(lattice_profile));
    lap->leg = leg;
    lap->Nsite = Nsite;
    lap->Nb = Nb;
    lap->bond2index = (int*)malloc(sizeof(int)*nl*Nb);
    lap->bondst = (double*)malloc(sizeof(double)*Nb);

    return lap;
}

void destroy_lattice_profile(lattice_profile* lap){
    free(lap->bond2index);
    free(lap->bondst);
    free(lap);
}

void lattice_profile_set_beta(lattice_profile* lap, double beta){
    lap->beta=beta;
}

system_state* create_system_state(int leg, int Nsite, int length){
    if(leg%2==1){
        printf("data_struct.c create_system_state : leg should be even number\n");
        exit(-1);
    }

    system_state* state = (system_state*)malloc(sizeof(system_state));
    state->leg = leg;
    state->Nsite = Nsite;
    state->length = length;
    state->sigma = (int*)malloc(sizeof(int)*Nsite);
    state->sequence = (int*)malloc(sizeof(int)*length);
    state->rng = gsl_rng_alloc(gsl_rng_mt19937);

    state->noo = 0;
    for(int i=0;i<length;++i) state->sequence[i]=-1;

    return state;
}

void destroy_system_state(system_state* state){
    free(state->sigma);
    free(state->sequence);
    gsl_rng_free(state->rng);
    free(state);
}

void system_state_set_rng(system_state* state, int seed){
    gsl_rng_set(state->rng,seed);
}

placeholder* create_placeholder(int leg, int Nsite, int length){
    if(leg%2==1){
        printf("data_struct.c create_placeholder : leg should be even number\n");
        exit(-1);
    }

    placeholder* ws = (placeholder*)malloc(sizeof(placeholder));
    ws->leg = leg;
    ws->Nsite = Nsite;
    ws->length = length;
    ws->sigma  = (int*)malloc(sizeof(int)*Nsite);
    ws->vfirst = (int*)malloc(sizeof(int)*Nsite);
    ws->vlast  = (int*)malloc(sizeof(int)*Nsite);
    ws->linkv  = (int*)malloc(sizeof(int)*length*leg);

    return ws;
}

void destroy_placeholder(placeholder* ws){
    free(ws->sigma);
    free(ws->vfirst);
    free(ws->vlast);
    free(ws->linkv);
    free(ws);
}

void adjust_cutoff(system_state* state, placeholder* ph, double buffer){
    if(state->noo*buffer>state->length){
        int newl = (int)state->noo*buffer;
        int leg    = ph->leg;
        int* newsq = (int*)malloc(sizeof(int)*newl);
        int* newlv = (int*)malloc(sizeof(int)*newl*leg);
        int length = state->length;
        for(int i=0;i<length;++i) newsq[i]=state->sequence[i];
        for(int i=length;i<newl;++i) newsq[i]=-1;

        free(state->sequence);
        state->sequence = newsq;
        state->length = newl;

        free(ph->linkv);
        ph->linkv = newlv;
        ph->length = newl;
    }
}

void sequence_doubling(lattice_profile* lap, system_state* state, placeholder* ws){
    if(state->leg!=ws->leg){
        printf("data_struct.c sequence_doubling : The number of leg in state and ws should be the same!\n");
        exit(-1);
    }
    else if(state->Nsite!=ws->Nsite){
        printf("data_struct.c sequence_doubling : The Nsite in state and ws should be the same!\n");
        exit(-1);
    }
    else if(state->length!=ws->length){
        printf("data_struct.c sequence_doubling : The length in state and ws should be the same!\n");
        exit(-1);
    }
    else if(state->leg!=lap->leg){
        printf("data_struct.c sequence_doubling : The number of leg in state and lap should be the same!\n");
        exit(-1);
    }
    int newl = state->length*2;
    int* newsq = (int*)malloc(sizeof(int)*newl);
    int length = state->length;
    int leg = state->leg;
    for(int i=0;i<length;++i){
        newsq[i]=state->sequence[i];
        newsq[i+state->length]=state->sequence[i];
    }

    free(state->sequence);
    state->sequence = newsq;
    state->length = newl;
    state->noo *= 2;

    free(ws->linkv);
    ws->linkv = (int*)malloc(sizeof(int)*newl*leg);
    ws->length = newl;

    lap->beta *=2;
}

#if 0
#define TEST_DATA_STRUCT
#endif
#ifdef TEST_DATA_STRUCT
int main(){
    int leg = 4;
    int Nsite = 32;
    int Nb = 64;
    int length = 1024;
    double beta = 20;

    for(int i=0;i<10000000;++i){
    lattice_profile* lap = create_lattice_profile(leg,Nsite,Nb);
    system_state* state = create_system_state(leg,Nsite,length);
    placeholder* ph = create_placeholder(leg,Nsite,length);

    lattice_profile_set_beta(lap,beta);
    sequence_doubling(lap,state,ph);

    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);
    }
    return 0;
}

estimators* create_estimators(int nsample, int nobs){
    estimators* est = (estimators*)malloc(sizeof(estimators));
    est->nsample=nsample;
    est->nobs=nobs;
    est->data=(double*)malloc(sizeof(double)*nobs*nsample);
    est->means=(double*)malloc(sizeof(double)*nobs);

    return est;
}

void destroy_estimators(estimators* est){
    free(est->data);
    free(est->means);
    free(est);
}

double estimators_get_data(const estimators* est, int i_obs, int i_sample){
    int nobs = est->nobs;
    return est->data[nobs*i_sample+i_obs];
}

void estimators_write_data(estimators* est, int i_obs, int i_sample, double data){
    int nobs = est->nobs;
    est->data[nobs*i_sample+i_obs]=data;
}

#endif
