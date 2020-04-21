#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "data_struct.h"
#include "estimator.h"
#include "heisenberg_model.h"

int MEASUREMENT_DILUTED_BILAYER = 0;
double* struct_factor_pi;
void measurement_diluted_bilayor(estimators* est, placeholder* ph, const system_state* state, const lattice_profile* lap, int nx, int ny, int i_sample){
    //int nobs    = est->nobs;
    int Nsite   = lap->Nsite;
    int length  = state->length;
    double beta = lap->beta;

    if(MEASUREMENT_DILUTED_BILAYER==0){
        MEASUREMENT_DILUTED_BILAYER = 1;
        struct_factor_pi = (double*)malloc(sizeof(double)*Nsite);
        for(int i=0;i<Nsite;++i){
            int z = i/(nx*ny);
            int x = (i%(nx*ny))%nx;
            int y = (i%(nx*ny))/nx;
            struct_factor_pi[i] = 1-((x+y+z)%2)*2;
        }
    }

    double mz=0,msz=0;
    int vol=0;
    int i,index;

    for(i=0;i<Nsite;++i){
        ph->sigma[i] = state->sigma[i];
        mz += state->sigma[i];
        msz += state->sigma[i]*struct_factor_pi[i];
        if(state->sigma[i]!=0) vol++;
    }
    double mag  = mz/Nsite*0.5;
    double mag2 = mz*mz/Nsite*0.25*beta;
    double noo  = state->noo;
    double noo2 = noo*noo;

    double ms1=0;
    double ms2=0;
    double ms4=0;
    int sp,i_bond,type;
    for(int p=0;p<length;++p){
        sp = state->sequence[p];
        i_bond = sp/2;
        type = sp%2;
        if(type==1){
            index = lap->bond2index[i_bond];
            msz += -4*(ph->sigma[index])*struct_factor_pi[index];
        }
        ms1 += fabs(msz);
        ms2 += msz*msz;
        ms4 += msz*msz*msz*msz;
        hm_propagate_state(ph,lap,state->sequence[p]);
    }

    ms1 = ms1*0.5/length/Nsite;
    ms2 = ms2*0.25/length/(Nsite*Nsite);
    ms4 = ms4*0.125/length/(Nsite*Nsite*Nsite*Nsite);

    estimators_write_data(est,0,i_sample,mag);
    estimators_write_data(est,1,i_sample,mag2);
    estimators_write_data(est,2,i_sample,noo);
    estimators_write_data(est,3,i_sample,noo2);
    estimators_write_data(est,4,i_sample,ms1);
    estimators_write_data(est,5,i_sample,ms2);
    estimators_write_data(est,6,i_sample,ms4);
    estimators_write_data(est,7,i_sample,vol);
}


#if 0
#define TEST_OBSERVABLES_DILUTED_BILAYER_HEISENBERG
#endif
#ifdef TEST_OBSERVABLES_DILUTED_BILAYER_HEISENBERG

#include "diluted_bilayer_heisenberg.h"
#include "heisenberg_model.h"

int main(){
    int seed=89798332;
    int nx=20;
    int ny=20;
    double p=0.407;
    double jbond=0.15;
    double beta=20;
    int nthermal=10000;
    int nsweep=20000;
    int nblock=20;

    lattice_profile* lap;
    system_state* state;
    placeholder* ph;
    estimators* est = create_estimators(nsweep,4);

    set_site_diluted_bilayer(&lap,&state,&ph,nx,ny,jbond,p,seed);

    lattice_profile_set_beta(lap,beta);

    for(int i=0;i<nthermal;++i){
        hm_monte_carlo_sweep(state,ph,lap);
        adjust_cutoff(state,ph,1.5);
    }

    for(int b=0;b<nblock;++b){
        for(int i=0;i<nsweep;++i){
            hm_monte_carlo_sweep(state,ph,lap);
            measurement_diluted_bilayor(est,ph,state,lap,nx,ny,i);
        }
        estimators_save_data(est,lap,"test.txt");
    }


    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);

    return 0;
}
#endif
