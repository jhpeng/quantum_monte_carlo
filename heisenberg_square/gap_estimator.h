#ifndef gap_estimator_h
#define gap_estimator_h


void gap_estimator_setup_workspace(int length);

void gap_estimator_free_memory();

void gap_estimator_collect_data(
                int* sequence, 
                int length, 
                int* sigma0, 
                int* sigmap, 
                int nsite, 
                int* bond2index, 
                double* structfactor);

void gap_estimator_calc_power_spectrum();

void gap_estimator_ave_power_spectrum_fileout(char* prefix);


#endif
