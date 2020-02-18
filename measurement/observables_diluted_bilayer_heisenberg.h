#ifndef observables_diluted_bilayer_heisenberg_h
#define observables_diluted_bilayer_heisenberg_h

#include "data_struct.h"

void measurement_diluted_bilayor(
            estimators* est, 
            placeholder* ph, 
            const system_state* state, 
            const lattice_profile* lap, 
            int nx, 
            int ny, 
            int i_sample);

#endif
