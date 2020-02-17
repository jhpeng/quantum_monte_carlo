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

static void hm_construct_link_vertex_list(placeholder* ph,const lattice_profile* lap, const system_state* state){
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
