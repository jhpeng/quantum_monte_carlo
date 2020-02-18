#include <stdio.h>

#include "data_struct.h"
#include "estimator.h"

void measurement_diluted_bilayor(estimators* est, placeholder* ph, const system_state* state, const lattice_profile* lap, int i_sample){
    int nobs  = est->nobs;
    int Nsite = lap->Nsite;

    double mz =0;

    for(int i=0;i<Nsite;++i){
        mz += state->sigma[i];
    }
    double mag  = mz/Nsite;
    double mag2 = mz*mz/Nsite;
    double noo  = state->noo;
    double noo2 = noo*noo;

    estimators_write_data(est,0,i_sample,mag);
    estimators_write_data(est,1,i_sample,mag2);
    estimators_write_data(est,2,i_sample,noo);
    estimators_write_data(est,3,i_sample,noo2);
}


#if 1
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

    for(int i=0;i<nsweep;++i){
        hm_monte_carlo_sweep(state,ph,lap);

        
        measurement_diluted_bilayor(est,ph,state,lap,i);
    }

    estimators_save_data(est,lap,"test.txt");

    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);

    return 0;
}
#endif
