#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/*
**  define global variables here
*/

/* The flip rule is the instruction for flipping operator.
**   | 0 | 1 | 2 | 3 | 4 | 5
------------------------------
** a | 1 | 0 | 4 | 5 | 2 | 3
** b | 0 | 1 | 3 | 2 | 5 | 4
** c | 1 | 0 | 4 | 5 | 2 | 3
** d | 0 | 1 | 3 | 2 | 5 | 4
**
**
**  |    |        |    |
** c|   c|       d|   d| 
** --------      --------
** |\\\\\\|------|\\\\\\|
** --------      --------
**  |    |        |    |
** a|   a|       b|   b|
**
*/
static int flip_rule[24]={
    1,0,4,5,2,3,
    0,1,3,2,5,4,
    1,0,4,5,2,3,
    0,1,3,2,5,4};

/*
** L   : total length of the operator sequence
** Noo : number of non-identity operator in the operator sequence
** Sequence : pointer to the operator sequence
** Linkv    : pointer to the link-vertex list (8L)
*/
int L,Noo;
int* Sequence;
int* Linkv;

/*
** Nsite : the total number of the spin
** Nc    : the number of spin in the cluster
** Sigma0 : the pointer to the current state
** Sigmap : the pointer to the propagate state
** Vfirst : the pointer to the first list
** Vlast  : the pointer to the last list
*/
int Nsite;
int Nc;
int* Sigma0;
int* Sigmap;
int* Vfirst;
int* Vlast;

/*
** Nj : the total number of the J-bonds
** Nq : the totla number of the Q-bond-pairs
** Bond2index : the pointer which mapping the bond on four spin site
** Bondst : the pointer to the bond strenght
*/
int Nx,Ny,Nz;
int Nj,Nq;
int* Bond2index;
double* Bondst;
double* StructFactor;

/*
** Nobs    : the total number of the observables
** Nsample : the total number of the Monte Carlo sample
** Data : the pointer to the collected data
*/
int Nobs,Nsample;
double* Data;
char Filename[128]="test.txt";

/*
** Nblock  : the total number of the block data
** Beta    : inverse temperature
** rng     : gsl_rng
** Mode    : the mode for calculate observable
**             0 -> normal scheme
**             1 -> beta doubling
*/
int Nblock,Nther,Seed;
double Beta,Jbond,Qbond,P;
gsl_rng* rng;
int Mode,LatticeType;
int Nit;


/* -------------------------------------------------- **
** ---------------- SSE algorithm ------------------- **
** -------------------------------------------------- */

void propagate_state(int sp){
    if(sp==-1) return;

    int type   = sp%6;
    if(type==0 || type==2) return;
    
    int i_bond = sp/6;
    if(type==1 || type==4){
        Sigmap[Bond2index[i_bond*4+0]] *= -1;
        Sigmap[Bond2index[i_bond*4+1]] *= -1;
    }
    else if(type==3){
        Sigmap[Bond2index[i_bond*4+2]] *= -1;
        Sigmap[Bond2index[i_bond*4+3]] *= -1;
    }
    else if(type==5){
        Sigmap[Bond2index[i_bond*4+0]] *= -1;
        Sigmap[Bond2index[i_bond*4+1]] *= -1;
        Sigmap[Bond2index[i_bond*4+2]] *= -1;
        Sigmap[Bond2index[i_bond*4+3]] *= -1;
    }
}

void diagonal_update(){
    int i_bond,s1,s2,s3,s4;
    double dis;

    for(int i=0;i<Nsite;++i){
        Sigmap[i] = Sigma0[i];
    }

    for(int p=0;p<L;++p){
        if(Sequence[p]==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*2*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                    Sequence[p] = i_bond*6;
                    Noo++;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*4*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                    Sequence[p] = i_bond*6+2;
                    Noo++;
                }
            }
        }
        else if(Sequence[p]%6==0){
            i_bond = Sequence[p]/6;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<2*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else if(Sequence[p]%6==2){
            i_bond = Sequence[p]/6;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<4*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else propagate_state(Sequence[p]);
    }
}

void diagonal_update_exchange(){
    int i_bond,s1,s2,s3,s4;
    double dis;

    for(int i=0;i<Nsite;++i){
        Sigmap[i] = Sigma0[i];
    }

    for(int p=0;p<L;++p){
        if(Sequence[p]%6==0){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]){
                    Sequence[p] = i_bond*6;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]*2<Bondst[i_bond]){
                    Sequence[p] = i_bond*6+2;
                }
            }
        }
        else if(Sequence[p]%6==2){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]*2){
                    Sequence[p] = i_bond*6;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]){
                    Sequence[p] = i_bond*6+2;
                }
            }
        }
        else propagate_state(Sequence[p]);
    }
}

void construct_link_vertex_list(){
    for(int i=0;i<(8*L);++i) Linkv[i]=-1;
    for(int i=0;i<Nsite;++i){
        Vfirst[i]=-1;
        Vlast[i] =-1;
    }

    int i_bond,index,nu0,nu1;
    for(int p=0;p<L;++p){
        if(Sequence[p]!=-1){
            i_bond = Sequence[p]/6;
            for(int i=0;i<4;++i){
                index = Bond2index[i_bond*4+i];
                if(index!=-1){
                    nu0 = 8*p+i;
                    nu1 = Vlast[index];
                    if(nu1!=-1){
                        Linkv[nu0] = nu1;
                        Linkv[nu1] = nu0;
                    }
                    else Vfirst[index] = nu0;

                    Vlast[index] = nu0+4;
                }
            }
        }
    }

    for(int i=0;i<Nsite;++i){
        if(Vlast[i]!=-1){
            Linkv[Vlast[i]] = Vfirst[i];
            Linkv[Vfirst[i]] = Vlast[i];
        }
    }
}

void loop_update(){
   int nu0,nup,nun,flip=-1;

    for(nu0=0;nu0<(L*8);nu0+=2){
        if(Linkv[nu0]>=0){
            nun = nu0;
            if(gsl_rng_uniform_pos(rng)<0.5) flip=-1;
            else flip=-2;
            while(Linkv[nun]>=0){
                Linkv[nun] = flip;
                nup = nun^1;
                nun = Linkv[nup];
                Linkv[nup]=flip;
            }
        }
    } 
}

void flip_bit_operator(){
    int nu,i_flip_rule,type,i_bond,index;
    
    for(nu=0;nu<(8*L);nu+=2){
        if(Linkv[nu]==-2){
            type = Sequence[nu/8]%6;
            i_bond = Sequence[nu/8]/6;
            i_flip_rule = ((nu%8)/2)*6+type;
            Sequence[nu/8] = i_bond*6+flip_rule[i_flip_rule];
        }
    }

    for(index=0;index<Nsite;++index){
        nu = Vlast[index];
        if(nu==-1){
            if(gsl_rng_uniform_pos(rng)<0.5){
                Sigma0[index] *= -1;
            }
        }
        else if(nu>=0){
            if(Linkv[nu]==-2){
                Sigma0[index] *= -1;
            }
        }
    }
}


void beta_doubling(){
    int length =2*L;
    int* seq = (int*)malloc(length*sizeof(int));
    for(int p=0;p<length;++p){
        if(p<L) seq[p] = Sequence[p];
        else seq[p] = Sequence[p-L];
    }
    free(Sequence);
    free(Linkv);

    Sequence = seq;
    Linkv = (int*)malloc(8*length*sizeof(int));

    L     = length;
    Beta *= 2;
    Noo  *= 2;
}


/* --------------------------------------------------------- **
** ---------------- Setting the model ---------------------- **
** --------------------------------------------------------- */


void set_random_number(int seed){
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
}

void construct_cluster(int index, int* cluster, int* cs){
    int nc=0;
    int ic=0;
    int check=1;
    int k0=index;
    int k1;

    cluster[nc] = k0;
    nc++;
    
    while(check){
        k0 = cluster[ic];
        for(int i_bond=0;i_bond<Nj;++i_bond){
            int key1=0;
            if(Bond2index[4*i_bond+0]==k0) {
                k1  = Bond2index[4*i_bond+1];
                key1 = 1;
            }
            else if(Bond2index[4*i_bond+1]==k0) {
                k1=Bond2index[4*i_bond+0];
                key1 = 1;
            }

            if(key1){
                int key2=1;
                for(int i=0;i<nc;++i){
                    if(cluster[i]==k1){
                        key2=0;
                        i=nc;
                    }
                }

                if(key2){
                    cluster[nc] = k1;
                    nc++;
                }
            }
        }

        ic++;
        if(ic==nc) check=0;
    }

    *cs = nc;
}

void set_lattice_bilayer_site_diluted(int nx, int ny, double jbond, double p){
    int i,j,t,q,index1,index2;
    Nz    = 2;
    Nsite = 2*nx*ny;
    Nj    = 5*nx*ny;
    Nq    = 0;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    int* epsilon = (int*)malloc(2*nx*ny*sizeof(int));
    for(i=0;i<nx*ny;++i){
        if(gsl_rng_uniform_pos(rng)<p) {
            epsilon[i]=0;
            epsilon[i+nx*ny]=0;
        }
        else {
            epsilon[i]=1;
            epsilon[i+nx*ny]=1;
        }
    }

    int r_bond=0;
    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==2){
            index1 = i+nx*j+nx*ny;
            index2 = ((i+1)%nx)+nx*j+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==3){
            index1 = i+nx*j+nx*ny;
            index2 = i+nx*((j+1)%ny)+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==4){
            index1 = i+nx*j;
            index2 = i+nx*j+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = jbond;

                r_bond++;
            }
        }
    }

    Nj = r_bond;
    Nc = Nsite;
    for(i=0;i<2*nx*ny;++i){
        if(epsilon[i]==0){
            Sigma0[i] = 0;
            Nc--;
        }
        else if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }

    free(epsilon);
}

void set_lattice_bilayer_site_diluted_largest_cluster(int nx, int ny, double jbond, double p){
    int i,j,t,q,index1,index2;
    Nz    = 2;
    Nsite = 2*nx*ny;
    Nj    = 5*nx*ny;
    Nq    = 0;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    int* epsilon = (int*)malloc(2*nx*ny*sizeof(int));
    for(i=0;i<nx*ny;++i){
        if(gsl_rng_uniform_pos(rng)<p) {
            epsilon[i]=0;
            epsilon[i+nx*ny]=0;
        }
        else {
            epsilon[i]=1;
            epsilon[i+nx*ny]=1;
        }
    }

    int r_bond=0;
    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==2){
            index1 = i+nx*j+nx*ny;
            index2 = ((i+1)%nx)+nx*j+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==3){
            index1 = i+nx*j+nx*ny;
            index2 = i+nx*((j+1)%ny)+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = 1;

                r_bond++;
            }
        }
        else if(q==4){
            index1 = i+nx*j;
            index2 = i+nx*j+nx*ny;
            if(epsilon[index1] && epsilon[index2]){
                Bond2index[r_bond*4+0] = index1;
                Bond2index[r_bond*4+1] = index2;
                Bond2index[r_bond*4+2] = -1;
                Bond2index[r_bond*4+3] = -1;
                Bondst[r_bond] = jbond;

                r_bond++;
            }
        }
    }

    Nj = r_bond;

    int* cluster = (int*)malloc(Nsite*sizeof(int));
    int cs=0;
    int max_s=0;
    int index;
    for(i=0;i<Nsite && max_s<nx*ny;++i){
        if(epsilon[i]){
            construct_cluster(i,cluster,&cs);
            if(cs>max_s){
                max_s=cs;
                index=i;
            }
        }
    }

    construct_cluster(index,cluster,&cs);

    for(i=0;i<Nsite;++i) Sigma0[i]=0;
    for(i=0;i<cs;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[cluster[i]]=1;
        else Sigma0[cluster[i]]=-1;
    }

    r_bond=0;
    for(i=0;i<Nj;++i){
        index1 = Bond2index[4*i+0]; 
        index2 = Bond2index[4*i+1];
        int key=0;
        for(index=0;index<cs;++index){
            if(index1==cluster[index]) key=1;
            else if(index2==cluster[index]) key=1;
        }

        if(key){
            Bond2index[4*r_bond+0] = Bond2index[4*i+0];
            Bond2index[4*r_bond+1] = Bond2index[4*i+1];
            Bondst[r_bond] = Bondst[i];
            r_bond++;
        }
    }

    Nj = r_bond;
    Nc = cs;

    free(epsilon);
    free(cluster);
}

void create_structure_factor(int nx, int ny, int nz, int nsite){
    int i,x,y,z;
    StructFactor = (double*)malloc(nsite*sizeof(double));
    
    for(i=0;i<nx*ny*nz;++i){
        x = (i%(nx*ny))%nx;
        y = (i%(nx*ny))/nx;
        z = i/(nx*ny);

        if((x+y+z)%2) StructFactor[i] = -1;
        else StructFactor[i] = 1;
    }
}

void set_sequence_length(int length){
    if(Sequence==NULL){
        Sequence = (int*)malloc(length*sizeof(int));
        Noo=0;
        for(int p=0;p<length;++p){
            Sequence[p]=-1;
        }
    }
    else{
        int* temp = (int*)malloc(length*sizeof(int));
        for(int p=0;p<length;++p){
            if(p<L) temp[p] = Sequence[p];
            else temp[p]=-1;
        }
        free(Sequence);
        Sequence = temp;
    }

    if(Linkv!=NULL) free(Linkv);
    Linkv = (int*)malloc(8*length*sizeof(int));

    L = length;
}

void set_estimator(int n_obs, int n_sample, int n_block){
    Data = (double*)malloc(n_obs*n_sample*sizeof(double));
    Nobs = n_obs;
    Nsample = n_sample;
    Nblock = n_block;
}

void free_memory(){
    free(Sequence);
    free(Linkv);
    free(Sigma0);
    free(Sigmap);
    free(Vfirst);
    free(Vlast);
    free(Bond2index);
    free(Bondst);
    free(StructFactor);
    free(Data);
}



/* --------------------------------------------------------- **
** ----------------------- Estimator ----------------------- **
** --------------------------------------------------------- */

void calc_mean(double* mean){
    for(int i_obs=0;i_obs<Nobs;++i_obs){
        mean[i_obs]=0;
        for(int i=0;i<Nsample;++i) mean[i_obs]+=Data[Nobs*i+i_obs];
        mean[i_obs] = mean[i_obs]/Nsample;
    }
}

void estimator_stdout(){
    double mean[Nobs];
    calc_mean(mean);

    printf("Nsample : %d\n",Nsample);
    printf("beta    : %.5e\n",Beta);
    for(int i=0;i<Nobs;++i) printf("%.5e ",mean[i]);
    printf("\n---------------------------------------\n");
}

void estimator_fileout(char* filename){
    double mean[Nobs];
    calc_mean(mean);

    FILE* ofile = fopen(filename,"a");
    fprintf(ofile,"%.5e ",Beta);
    for(int i=0;i<Nobs;++i) fprintf(ofile,"%.5e ",mean[i]);
    fprintf(ofile,"\n");
    fclose(ofile);
}

void measure_with_propagate_state(int i_sample){
    double m1,m2,mz=0;
    double ms=0;
    double ms1=0;
    double ms2=0;
    double ms4=0;
    for(int i=0;i<Nsite;++i) Sigmap[i]=Sigma0[i];
    for(int i=0;i<Nsite;++i) mz+=Sigma0[i];
    for(int i=0;i<Nsite;++i) ms+=Sigma0[i]*StructFactor[i];

    m1 = mz*0.5;
    m2 = mz*mz*0.25;
    Data[Nobs*i_sample+0] = m1/Nc;
    Data[Nobs*i_sample+1] = m2*Beta/Nc;

    for(int p=0;p<L;++p){
        int sp = Sequence[p];
        int type = sp%6;
        int i_bond = sp/6;
        int i = Bond2index[i_bond*4+0];

        if(type==1) ms += -4*Sigmap[i]*StructFactor[i];

        ms1 += fabs(ms);
        ms2 += ms*ms;
        ms4 += ms*ms*ms*ms;
        propagate_state(sp);
    }

/*-------------------- check inner product ---------------------*/
    int check = 0;
    for(int i=0;i<Nsite;++i){
        if(Sigma0[i]!=Sigmap[i]) check=1;
    }
    if(check){
        printf("inner product become zero!\n");
        exit(-1);
    }
/*--------------------------------------------------------------*/

    ms1 = ms1/L/Nc*0.5;
    ms2 = ms2/L/Nc/Nc*0.25;
    ms4 = ms4/L/Nc/Nc/Nc/Nc*0.0625;

    Data[Nobs*i_sample+2] = ms1;
    Data[Nobs*i_sample+3] = ms2;
    Data[Nobs*i_sample+4] = ms4;

    double noo = Noo;
    Data[Nobs*i_sample+5] = noo;
    Data[Nobs*i_sample+6] = noo*noo;

    double se = 0;
    for(int i=0;i<Nj;++i) se += Bondst[i];
    Data[Nobs*i_sample+7] = (-noo/Beta + se*0.25)/Nc;
}

/* ----------------------------------------------- **
** ------------------ getopt --------------------- **
** ----------------------------------------------- */ 

int Help;
void set_opt(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hx:y:j:b:n:k:t:e:s:f:m:l:p:"))!=-1){
        switch(c){
            case 'h':
                Help=1;
                printf("usage:\n");
                printf("\t-h print this help\n");
                printf("\t-l lattice type for the simulation\n");
                printf("\t\t 0 : 2d ladder\n");
                printf("\t\t 1 : 2d herringbone\n");
                printf("\t-m mode for calculate observable\n");
                printf("\t\t 0 : normal scheme\n");
                printf("\t\t 1 : beta-doubling scheme\n");
                printf("\t\t 2 : beta increasing scheme\n");
                printf("\t\t 3 : beta increasing scheme without propagate state\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-j <Jp/J ratio> default 1.0\n");
                printf("\t-p <strength of randomness> default 0.0\n");
                printf("\t-b <beta> default 4.0\n");
                printf("\t-n <Nsample> default 2000\n");
                printf("\t-k <Nblock> default 50\n");
                printf("\t-t <Nther> default 2000\n");
                printf("\t-e number of iteration for beta-doubing scheme\n");
                printf("\t-s <seed of random number generator> default 1\n");
                printf("\t-f <the file name of output data> default \"test.txt\"\n");
                break;
            case 'l':
                LatticeType=atoi(optarg);
                break;
            case 'm':
                Mode=atoi(optarg);
                break;
            case 'x':
                Nx=atoi(optarg);
                break;
            case 'y':
                Ny=atoi(optarg);
                break;
            case 'j':
                Jbond=atof(optarg);
                break;
            case 'p':
                P=atof(optarg);
                break;
            case 'b':
                Beta=atof(optarg);
                break;
            case 'n':
                Nsample=atoi(optarg);
                break;
            case 'k':
                Nblock=atoi(optarg);
                break;
            case 't':
                Nther=atoi(optarg);
                break;
            case 'e':
                Nit=atoi(optarg);
                break;
            case 's':
                Seed=atoi(optarg);
                break;
            case 'f':
                strcpy(Filename,optarg);
                break;
        }
    }
}


/* ----------------------------------------------- **
** ------------------- main ---------------------- **
** ----------------------------------------------- */ 

int main(int argc, char** argv){
    int length=100;
    int n_obs=4;
    double buffer=1.3;

    /*--------------default value----------------*/
    Beta = 4096;
    Seed = 9237912;
    Nx = 48;
    Ny = 48;
    Nz = 1;
    Jbond = 1.0;
    P = 0;
    Nther = 20000;
    Nsample = 2000;
    Nblock = 50;
    Mode = 0;
    LatticeType = 0;
    Nit = 5;


    /*----------------get option-----------------*/
    set_opt(argc,argv);
    if(Help) return 0;


    /*-----------set random generator------------*/
    set_random_number(Seed);


    /*---------------set lattice-----------------*/
    if(LatticeType==0) set_lattice_bilayer_site_diluted(Nx,Ny,Jbond,P);
    else if(LatticeType==1) set_lattice_bilayer_site_diluted_largest_cluster(Nx,Ny,Jbond,P);
    create_structure_factor(Nx,Ny,Nz,Nsite);

    for(int i=0;i<Nj;++i){
        printf("%d %d \n",Bond2index[4*i+0],Bond2index[4*i+1]);
    }
    printf("Nsite = %d\n",Nsite);
    printf("Nc    = %d\n",Nc);
    printf("Nj    = %d\n",Nj);

    set_sequence_length(length);

    if(Mode==0 || Mode==1 || Mode==2) n_obs=8;
    else if(Mode==3) n_obs=4;
    set_estimator(n_obs,Nsample,Nblock);



/**************************************************************/
/*********************** Normal Scheme ************************/
/**************************************************************/
    if(Mode==0){
        /*---------------Thermalization--------------*/
        for(int i_sample=0;i_sample<Nther;++i_sample){
            diagonal_update();
            construct_link_vertex_list();
            loop_update();
            flip_bit_operator();
            if(Noo*buffer>L){
                length = (int)(Noo*buffer)+10;
                set_sequence_length(length);
            }
        }

        /*---------------Measurement-----------------*/
        for(int k=0;k<Nblock;++k){
            for(int i_sample=0;i_sample<Nsample;++i_sample){
                diagonal_update();
                construct_link_vertex_list();
                loop_update();
                flip_bit_operator();

                measure_with_propagate_state(i_sample);
            }
            estimator_fileout(Filename);
        }
    }
/**************************************************************/
/******************** Beta-doubling Scheme ********************/
/**************************************************************/
    else if(Mode==1){
        for(int it=0;it<2*Nit;++it){
            /*---------------Thermalization--------------*/
            for(int i_sample=0;i_sample<Nther;++i_sample){
                diagonal_update();
                construct_link_vertex_list();
                loop_update();
                flip_bit_operator();
                if(Noo*buffer>L){
                    length = (int)(Noo*buffer)+10;
                    set_sequence_length(length);
                }
            }

            /*---------------Measurement-----------------*/
            for(int k=0;k<Nblock;++k){
                for(int i_sample=0;i_sample<Nsample;++i_sample){
                    diagonal_update();
                    construct_link_vertex_list();
                    loop_update();
                    flip_bit_operator();

                    measure_with_propagate_state(i_sample);
                }
                estimator_fileout(Filename);
            }
            if(it%2==1){
                beta_doubling();
            }
        }
    }

    free_memory();
    return 0;
}
