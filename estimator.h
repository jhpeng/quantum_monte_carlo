#ifndef estimator_h
#define estimator_h

#include <stdio.h>
#include "data_struct.h"

void estimator_save_data(
            estimators* est, 
            const lattice_profile* lap, 
            const char* filename);

#endif
