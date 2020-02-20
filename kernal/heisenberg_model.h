#ifndef heisenberg_model_h
#define heisenberg_model_h

#include "data_struct.h"

void hm_propagate_state(
            placeholder* ph, 
            const lattice_profile* lap, 
            int sp);

void hm_monte_carlo_sweep(
            system_state* state, 
            placeholder* ph, 
            const lattice_profile* lap);

#endif
