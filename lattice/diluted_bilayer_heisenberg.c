#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#include "data_struct.h"


void set_site_diluted_bilayer(lattice_profile** lap, system_state** state, placeholder** ph, int nx, int ny, double jbond, double p, int seed){
    int i,j,t,q,index1,index2;

    int leg    = 4;
    int Nsite  = 2*nx*ny;
    int Nb     = 5*nx*ny;
    int length = 1024;

    *lap = create_lattice_profile(leg,Nsite,Nb);
    *state = create_system_state(leg,Nsite,length);
    *ph = create_placeholder(leg,Nsite,length);

    system_state_set_rng(*state,seed);

    int* epsilon = (int*)malloc(Nsite*sizeof(int));
    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos((*state)->rng)<p) epsilon[i] = 0;
        else epsilon[i] = 1;
    }

    for(int i_bond=0;i_bond<Nb;++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]) (*lap)->bondst[i_bond] = 1;
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]) (*lap)->bondst[i_bond] = 1;
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==2){
            index1 = i+nx*j+nx*ny;
            index2 = ((i+1)%nx)+nx*j+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]) (*lap)->bondst[i_bond] = 1;
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==3){
            index1 = i+nx*j+nx*ny;
            index2 = i+nx*((j+1)%ny)+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]) (*lap)->bondst[i_bond] = 1;
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==4){
            index1 = i+nx*j;
            index2 = i+nx*j+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]) (*lap)->bondst[i_bond] = jbond;
            else (*lap)->bondst[i_bond] = 0;
        }
    }
    for(i=0;i<Nsite;++i){
        if(epsilon[i]==0) (*state)->sigma[i]=0;
        else if(gsl_rng_uniform_pos((*state)->rng)<0.5) (*state)->sigma[i]=1;
        else (*state)->sigma[i]=-1;
    }

    free(epsilon);
}

void set_dimer_diluted_bilayer(lattice_profile** lap, system_state** state, placeholder** ph, int nx, int ny, double jbond, double p, int seed){
    int i,j,t,q,index1,index2;

    int leg    = 4;
    int Nsite  = 2*nx*ny;
    int Nb     = 5*nx*ny;
    int length = 1024;

    *lap = create_lattice_profile(leg,Nsite,Nb);
    *state = create_system_state(leg,Nsite,length);
    *ph = create_placeholder(leg,Nsite,length);

    system_state_set_rng(*state,seed);

    int* epsilon = (int*)malloc(Nsite*sizeof(int));
    for(i=0;i<nx*ny;++i){
        if(gsl_rng_uniform_pos((*state)->rng)<p){
            epsilon[i] = 0;
            epsilon[i+nx*ny] = 0;
        }
        else{
            epsilon[i] = 1;
            epsilon[i+nx*ny] = 1;
        }
    }

    int nb=0;
    int i_bond;
    int* rlist = malloc(sizeof(int)*Nb);
    for(i_bond=0;i_bond<Nb;++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]){
                (*lap)->bondst[i_bond] = 1;
                rlist[nb] = i_bond;
                nb++;
            }
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]){
                (*lap)->bondst[i_bond] = 1;
                rlist[nb] = i_bond;
                nb++;
            }
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==2){
            index1 = i+nx*j+nx*ny;
            index2 = ((i+1)%nx)+nx*j+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]){
                (*lap)->bondst[i_bond] = 1;
                rlist[nb] = i_bond;
                nb++;
            }
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==3){
            index1 = i+nx*j+nx*ny;
            index2 = i+nx*((j+1)%ny)+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]){
                (*lap)->bondst[i_bond] = 1;
                rlist[nb] = i_bond;
                nb++;
            }
            else (*lap)->bondst[i_bond] = 0;
        }
        else if(q==4){
            index1 = i+nx*j;
            index2 = i+nx*j+nx*ny;
            (*lap)->bond2index[i_bond*2+0] = index1;
            (*lap)->bond2index[i_bond*2+1] = index2;
            if(epsilon[index1] && epsilon[index2]){
                (*lap)->bondst[i_bond] = jbond;
                rlist[nb] = i_bond;
                nb++;
            }
            else (*lap)->bondst[i_bond] = 0;
        }
    }

    for(i=0;i<Nsite;++i){
        if(epsilon[i]==0) (*state)->sigma[i]=0;
        else if(gsl_rng_uniform_pos((*state)->rng)<0.5) (*state)->sigma[i]=1;
        else (*state)->sigma[i]=-1;
    }

    lattice_profile* new_lap = create_lattice_profile(leg,Nsite,nb);
    for(i=0;i<nb;++i){
        i_bond = rlist[i];
        new_lap->bond2index[i*2+0] = (*lap)->bond2index[i_bond*2+0];
        new_lap->bond2index[i*2+1] = (*lap)->bond2index[i_bond*2+1];
        new_lap->bondst[i] = (*lap)->bondst[i_bond];
    }

    destroy_lattice_profile(*lap);
    *lap = new_lap;

    free(epsilon);
    free(rlist);
}


#if 0
#define TEST_DILUTED_BILAYER_HEISENBERG
#endif
#ifdef TEST_DILUTED_BILAYER_HEISENBERG
int main(){
    int seed=89798332;
    int nx=4;
    int ny=4;
    double p=0.1;
    double jbond=0.3;
    double beta=20;

    lattice_profile* lap;
    system_state* state;
    placeholder* ph;

    set_site_diluted_bilayer(&lap,&state,&ph,nx,ny,jbond,p,seed);

    lattice_profile_set_beta(lap,beta);

    for(int i=0;i<lap->Nb;++i)
        printf("bonds : %d %d %d %f\n",i,lap->bond2index[2*i+0],lap->bond2index[2*i+1],lap->bondst[i]);

    for(int i=0;i<lap->Nsite;++i)
        printf("sigma : %d %d\n",i,state->sigma[i]);

    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);

    return 0;
}
#endif
