#ifndef momentum_h
#define momentum_h

void momentum_setup_workspace(
                int nx, 
                int ny, 
                int nz, 
                int nk, 
                double* wave_vector, 
                int n_sample);

void momentum_free_memory();

void momentum_collect_data(
                int* sequence, 
                int length, 
                int* sigma0, 
                int* sigmap, 
                int nsite, 
                int* bond2index, 
                double beta, 
                int i_sample);

void momentum_calc_mean_fileout(char* prefix);

#endif
