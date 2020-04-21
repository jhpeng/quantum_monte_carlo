#ifndef scheme_diluted_bilayer_heisenberg_h
#define scheme_diluted_bilayer_heisenberg_h

void fine_temp_diluted_bilayer_heisenberg(
    int nx, int ny, 
    double p, 
    double jbond, 
    double beta_i, 
    double beta_f, 
    double beta_v, 
    int nther, 
    int nsample, 
    int seed, 
    const char* filename);
#endif
