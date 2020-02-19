#include "data_struct.h"
#include "estimator.h"
#include "diluted_bilayer_heisenberg.h"
#include "heisenberg_model.h"
#include "observables_diluted_bilayer_heisenberg.h"


void fine_temp_diluted_bilayer_heisenberg(int nx, int ny, double p, double jbond, double beta_i, double beta_f, double beta_v, int nther, int nsample, int seed, const char* filename){
    lattice_profile* lap;
    system_state* state;
    placeholder* ph;
    estimators* est = create_estimators(nsample,5);

    //set_site_diluted_bilayer(&lap,&state,&ph,nx,ny,jbond,p,seed);
    set_dimer_diluted_bilayer(&lap,&state,&ph,nx,ny,jbond,p,seed);

    double beta;
    for(beta=beta_i;beta<beta_f;beta+=beta_v){
        lattice_profile_set_beta(lap,beta);

        for(int i=0;i<nther;++i){
            hm_monte_carlo_sweep(state,ph,lap);
            adjust_cutoff(state,ph,1.5);
        }

        for(int i=0;i<nsample;++i){
            hm_monte_carlo_sweep(state,ph,lap);
            measurement_diluted_bilayor(est,ph,state,lap,nx,ny,i);
        }
        estimators_save_data(est,lap,filename);
    }

    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);
    destroy_estimators(est);
}
