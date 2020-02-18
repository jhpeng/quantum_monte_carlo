#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#include "data_struct.h"

void hm_propagate_state(placeholder* ph, const lattice_profile* lap, int sp){
    if(sp==-1) return;

    int type = sp%2;
    if(type==0) return;
    else{
        int i_bond = sp/2;
        ph->sigma[lap->bond2index[i_bond*2+0]] *= -1;
        ph->sigma[lap->bond2index[i_bond*2+1]] *= -1;
    }
}

static void hm_check_weight(placeholder* ph, const lattice_profile* lap, const system_state* state){
    int i;
    int Nsite   = lap->Nsite;
    int length  = state->length;
    
    for(i=0;i<Nsite;++i) ph->sigma[i] = state->sigma[i];
    for(i=0;i<length;++i) hm_propagate_state(ph,lap,state->sequence[i]);

    int check=1;
    for(i=0;i<Nsite && check;++i){
        if(ph->sigma[i]!=state->sigma[i]) check=0;
    }

    if(check==0){
        printf("something wrong!\n");
        exit(-1);
    }
}

static void hm_diagonal_update(system_state* state, placeholder* ph, const lattice_profile* lap){
    int i,p,i_bond,s1,s2;
    double dis;

    int Nsite   = lap->Nsite;
    int Nb      = lap->Nb;
    double beta = lap->beta;
    int length  = state->length;
    int noo     = state->noo;

    for(i=0;i<Nsite;++i) ph->sigma[i] = state->sigma[i];

    for(p=0;p<length;++p){
        if(state->sequence[p]==-1){
            i_bond = (int)(gsl_rng_uniform_pos(state->rng)*Nb);
            s1 = ph->sigma[lap->bond2index[i_bond*2+0]];
            s2 = ph->sigma[lap->bond2index[i_bond*2+1]];
            if(s1!=s2 && s1!=0){
                dis = gsl_rng_uniform_pos(state->rng);
                if(dis*2*(length-noo)<beta*lap->bondst[i_bond]*Nb){
                    state->sequence[p] = i_bond*2;
                    noo++;
                }
            }
        }
        else if(state->sequence[p]%2==0){
            i_bond = state->sequence[p]/2;
            dis = gsl_rng_uniform_pos(state->rng);
            if(beta*Nb*lap->bondst[i_bond]*dis<2*(length-noo+1)){
                state->sequence[p]=-1;
                noo--;
            }
        }
        else hm_propagate_state(ph,lap,state->sequence[p]);
    }

    state->noo = noo;
}

static void hm_construct_link_vertex_list(placeholder* ph, const lattice_profile* lap, const system_state* state){
    int Nsite   = lap->Nsite;
    int length  = state->length;

    for(int i=0;i<(4*length);++i) ph->linkv[i]=-1;
    for(int i=0;i<Nsite;++i){
        ph->vfirst[i]=-1;
        ph->vlast[i]=-1;
    }

    int i,p,i_bond,index,nu0,nu1;
    for(p=0;p<length;++p){
        if(state->sequence[p]!=-1){
            i_bond = state->sequence[p]/2;
            for(i=0;i<2;++i){
                index = lap->bond2index[i_bond*2+i];
                nu0 = 4*p+i;
                nu1 = ph->vlast[index];
                if(nu1!=-1){
                    ph->linkv[nu0]=nu1;
                    ph->linkv[nu1]=nu0;
                }
                else ph->vfirst[index]=nu0;

                ph->vlast[index]=nu0+2;
            }
        }
    }

    for(i=0;i<Nsite;++i){
        if(ph->vlast[i]!=-1){
            ph->linkv[ph->vlast[i]] = ph->vfirst[i];
            ph->linkv[ph->vfirst[i]] = ph->vlast[i];
        }
    }
}

static void hm_loop_update(placeholder* ph, system_state* state){
    int nu0,nup,nun,flip=-1;
    int length = ph->length;

    for(nu0=0;nu0<(length*4);nu0+=2){
        if(ph->linkv[nu0]>=0){
            nun = nu0;
            if(gsl_rng_uniform_pos(state->rng)<0.5) flip=-1;
            else flip=-2;
            while(ph->linkv[nun]>=0){
                ph->linkv[nun] = flip;
                nup = nun^1;
                nun = ph->linkv[nup];
                ph->linkv[nup]=flip;
            }
        }
    } 
}

static void hm_flip_bit_operator(placeholder* ph, system_state* state){
    int nu,type,i_bond,index;
    int length = ph->length;
    int Nsite  = ph->Nsite;
    
    for(nu=0;nu<(4*length);nu+=2){
        if(ph->linkv[nu]==-2){
            type = state->sequence[nu/4]%2;
            i_bond = state->sequence[nu/4]/2;
            state->sequence[nu/4] = i_bond*2+type^1;
        }
    }

    for(index=0;index<Nsite;++index){
        nu = ph->vlast[index];
        if(nu==-1){
            if(gsl_rng_uniform_pos(state->rng)<0.5) state->sigma[index]*= -1;
        }
        else if(ph->linkv[nu]==-2) state->sigma[index]*= -1;
    }
}

void hm_monte_carlo_sweep(system_state* state, placeholder* ph, const lattice_profile* lap){
    hm_diagonal_update(state,ph,lap);
    hm_construct_link_vertex_list(ph,lap,state);
    hm_loop_update(ph,state);
    hm_flip_bit_operator(ph,state);
    //hm_check_weight(ph,lap,state);
}

#if 0
#define TEST_HEISENBERG_MODEL
#endif
#ifdef TEST_HEISENBERG_MODEL

#include "diluted_bilayer_heisenberg.h"

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

    set_site_diluted_bilayer(&lap,&state,&ph,nx,ny,jbond,p,seed);

    lattice_profile_set_beta(lap,beta);

    for(int i=0;i<nthermal;++i){
        hm_monte_carlo_sweep(state,ph,lap);
        adjust_cutoff(state,ph,1.5);
    }

    for(int i=0;i<nsweep;++i){
        hm_monte_carlo_sweep(state,ph,lap);

        double noo = state->noo;
        double noo2 = noo*noo;
        printf("sweep %d : noo %f | noo^2 %f\n",i,noo,noo2);

    }

    destroy_placeholder(ph);
    destroy_system_state(state);
    destroy_lattice_profile(lap);

    return 0;
}
#endif
