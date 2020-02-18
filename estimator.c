#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_struct.h"

static void calc_means(estimators* est){
    int nobs    = est->nobs;
    int nsample = est->nsample;

    for(int i_obs=0;i_obs<nobs;++i_obs){
        est->means[i_obs] = 0;
        for(int i_sample=0;i_sample<nsample;++i_sample){
            est->means[i_obs] += estimators_get_data(est,i_obs,i_sample);
        }
        est->means[i_obs] = est->means[i_obs]/nsample;
    }
}

void estimator_save_data(estimators* est, const lattice_profile* lap, const char* filename){
    calc_means(est);

    FILE* outfile = fopen(filename,"a");
    int nobs = est->nobs;
    fprintf(outfile,"%lf ",lap->beta);
    for(int i=0;i<nobs;++i){
        fprintf(outfile,"%.18e ",est->means[i]);
    }
    fprintf(outfile,"\n");
    fclose(outfile);
}
