#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "gap_estimator.h"
#include "momentum.h"

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
static int L,Noo;
static int* Sequence;
static int* Linkv;

/*
** Nsite : the total number of the spin
** Sigma0 : the pointer to the current state
** Sigmap : the pointer to the propagate state
** Vfirst : the pointer to the first list
** Vlast  : the pointer to the last list
*/
static int Nsite;
static int* Sigma0;
static int* Sigmap;
static int* Vfirst;
static int* Vlast;

/*
** Nj : the total number of the J-bonds
** Nq : the totla number of the Q-bond-pairs
** Bond2index : the pointer which mapping the bond on four spin site
** Bondst : the pointer to the bond strenght
*/
static int Nx,Ny,Nz;
static int Nj,Nq;
static int* Bond2index;
static double* Bondst;
static double* StructFactor;

/*
** Nobs    : the total number of the observables
** Nsample : the total number of the Monte Carlo sample
** Data : the pointer to the collected data
*/
static int Nobs,Nsample;
static double* Data;
static char Filename[128]="test.txt";

/*
** Nblock  : the total number of the block data
** Beta    : inverse temperature
** rng     : gsl_rng
** Mode    : the mode for calculate observable
**             0 -> normal scheme
**             1 -> beta doubling
*/
static int Nblock,Nther,Seed;
static double Beta,Jbond,Qbond,P;
static gsl_rng* rng;
static int Mode,LatticeType;
static int Nit;

static double PI = 3.141592653589793;


/* -------------------------------------------------- **
** ---------------- SSE algorithm ------------------- **
** -------------------------------------------------- */

void propagate_state(int* sigmap, const int* bond2index, int sp){
    if(sp==-1) return;

    int type   = sp%6;
    if(type==0 || type==2) return;
    
    int i_bond = sp/6;
    if(type==1 || type==4){
        sigmap[bond2index[i_bond*4+0]] *= -1;
        sigmap[bond2index[i_bond*4+1]] *= -1;
    }
    else if(type==3){
        sigmap[bond2index[i_bond*4+2]] *= -1;
        sigmap[bond2index[i_bond*4+3]] *= -1;
    }
    else if(type==5){
        sigmap[bond2index[i_bond*4+0]] *= -1;
        sigmap[bond2index[i_bond*4+1]] *= -1;
        sigmap[bond2index[i_bond*4+2]] *= -1;
        sigmap[bond2index[i_bond*4+3]] *= -1;
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
        else propagate_state(Sigmap,Bond2index,Sequence[p]);
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
        else seq[p] = Sequence[2*L-1-p];
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

void set_lattice_ladder(int nx, int ny, double jbond, double p, int x_open, int y_open){
    int i,j,t,q,index1,index2;
    Nz    = 1;
    Nsite = nx*ny;
    Nj    = 2*nx*ny;
    Nq    = 0;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            if(i%2==0){
                if(gsl_rng_uniform_pos(rng)<0.5) Bondst[i_bond] = 1+(jbond-1)*(1+p);
                else Bondst[i_bond] = 1+(jbond-1)*(1-p);
            }
            else Bondst[i_bond] = 1;

            if(x_open){
                if(i%nx==(nx-1)) Bondst[i_bond] = 0;
            }
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = 1;

            if(y_open){
                if(j%ny==(ny-1)) Bondst[i_bond] = 0;
            }
        }
    }
    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }
}

void set_lattice_herringbone(int nx, int ny, double jbond, double p){
    int i,j,t,q,index1,index2;
    double g;
    Nz    = 1;
    Nsite = nx*ny;
    Nj    = 2*nx*ny;
    Nq    = 0;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;
            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;

            if(gsl_rng_uniform_pos(rng)<0.5) g = 1+(jbond-1)*(1+p);
            else g = 1+(jbond-1)*(1-p);

            i = i%4; j = j%4;

            if(i==1 && j==0) Bondst[i_bond] = g;
            else if(i==2 && j==1) Bondst[i_bond] = g;
            else if(i==3 && j==2) Bondst[i_bond] = g;
            else if(i==0 && j==3) Bondst[i_bond] = g;
            else Bondst[i_bond] = 1;
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);
            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;

            if(gsl_rng_uniform_pos(rng)<0.5) g = 1+(jbond-1)*(1+p);
            else g = 1+(jbond-1)*(1-p);

            i = i%4; j = j%4;

            if(i==0 && j==0) Bondst[i_bond] = g;
            else if(i==1 && j==1) Bondst[i_bond] = g;
            else if(i==2 && j==2) Bondst[i_bond] = g;
            else if(i==3 && j==3) Bondst[i_bond] = g;
            else Bondst[i_bond] = 1;
        }
    }
    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }
}

void set_lattice_bilayer(int nx, int ny, double jbond){
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

    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%(nx*ny);
        q = i_bond/(nx*ny);
        i = t%nx;
        j = t/nx;

        if(q==0){
            index1 = i+nx*j;
            index2 = ((i+1)%nx)+nx*j;

            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = 1;
        }
        else if(q==1){
            index1 = i+nx*j;
            index2 = i+nx*((j+1)%ny);

            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = 1;
        }
        else if(q==2){
            index1 = i+nx*j+nx*ny;
            index2 = ((i+1)%nx)+nx*j+nx*ny;
            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = 1;
        }
        else if(q==3){
            index1 = i+nx*j+nx*ny;
            index2 = i+nx*((j+1)%ny)+nx*ny;

            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = 1;
        }
        else if(q==4){
            index1 = i+nx*j;
            index2 = i+nx*j+nx*ny;

            Bond2index[i_bond*4+0] = index1;
            Bond2index[i_bond*4+1] = index2;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = jbond;
        }
    }

    for(i=0;i<2*nx*ny;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }
}

void create_structure_factor(int nx, int ny, int nz, int nsite){
    int i,x,y,z;
    StructFactor = (double*)malloc(nsite*sizeof(double));
    
    for(i=0;i<nsite;++i){
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
    for(int i=0;i<Nobs;++i) fprintf(ofile,"%.16e ",mean[i]);
    fprintf(ofile,"\n");
    fclose(ofile);
}

void measure_with_propagate_state(int i_sample){
    int i,j,p,i_bond,type,sp;
    double m1,m2,mz=0;
    double ms=0;
    double msx=0;
    double ms1=0;
    double ms2=0;
    double ms4=0;
    for(i=0;i<Nsite;++i) Sigmap[i]=Sigma0[i];
    for(i=0;i<Nsite;++i) mz+=Sigma0[i];
    for(i=0;i<Nsite;++i) ms+=Sigma0[i]*StructFactor[i];

    m1 = mz*0.5;
    m2 = mz*mz*0.25;

    for(p=0;p<L;++p){
        sp = Sequence[p];

        if(sp!=-1){
            type = sp%6;
            i_bond = sp/6;
            i = Bond2index[i_bond*4+0];
            j = Bond2index[i_bond*4+1];

            if(type==1){
                ms += -2*Sigmap[i]*StructFactor[i];
                ms += -2*Sigmap[j]*StructFactor[j];
            }

            msx += ms;
            ms1 += fabs(ms);
            ms2 += ms*ms;
            ms4 += ms*ms*ms*ms;
            propagate_state(Sigmap,Bond2index,sp);
        }
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

    double noo = Noo;

    if(Noo!=0){
        msx = Beta*(msx*msx+ms2)/noo/(noo+1)/Nsite/Nsite*0.25;
        ms1 = ms1/noo/Nsite*0.5;
        ms2 = ms2/noo/Nsite/Nsite*0.25;
        ms4 = ms4/noo/Nsite/Nsite/Nsite/Nsite*0.0625;
    }
    else{
        msx=0;
        ms1 = fabs(ms)/Nsite*0.5;
        ms2 = ms*ms/Nsite/Nsite*0.25;
        ms4 = ms*ms*ms*ms/Nsite/Nsite/Nsite/Nsite*0.0625;
    }

    Data[Nobs*i_sample+0] = m1/Nsite;
    Data[Nobs*i_sample+1] = m2*Beta/Nsite;

    Data[Nobs*i_sample+2] = ms1;
    Data[Nobs*i_sample+3] = ms2;
    Data[Nobs*i_sample+4] = ms4;
    Data[Nobs*i_sample+5] = msx;

    Data[Nobs*i_sample+6] = noo;
    Data[Nobs*i_sample+7] = noo*noo;

    double se = 0;
    for(int i=0;i<Nj;++i) se += Bondst[i];
    Data[Nobs*i_sample+8] = (-noo/Beta + se*0.25)/Nsite;
}

/* ----------------------------------------------- **
** ------------------ getopt --------------------- **
** ----------------------------------------------- */ 

int Help;
void set_opt(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hx:y:j:b:n:k:t:e:s:o:m:l:p:"))!=-1){
        switch(c){
            case 'h':
                Help=1;
                printf("usage:\n");
                printf("\t-h print this help\n");
                printf("\t-l lattice type for the simulation\n");
                printf("\t\t 0 : 2d ladder\n");
                printf("\t\t 1 : 2d herringbone\n");
                printf("\t\t 2 : 2d bilayer\n");
                printf("\t-m mode for calculate observable\n");
                printf("\t\t 0 : <normal>\n");
                printf("\t\t 1 : <beta-doubling>\n");
                printf("\t\t 2 : <beta increasing>\n");
                printf("\t\t 3 : <beta-doubling> <gap estimator>\n");
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
                printf("\t-o <the file name of output data> default \"test.txt\"\n");
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
            case 'o':
                strcpy(Filename,optarg);
                break;
        }
    }
}


/* ----------------------------------------------- **
** ------------------- main ---------------------- **
** ----------------------------------------------- */ 


#if 1
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
    if(LatticeType==0) set_lattice_ladder(Nx,Ny,Jbond,P,0,1);
    else if(LatticeType==1) set_lattice_herringbone(Nx,Ny,Jbond,P);
    else if(LatticeType==2) set_lattice_bilayer(Nx,Ny,Jbond);
    create_structure_factor(Nx,Ny,Nz,Nsite);

    set_sequence_length(length);

    n_obs=9;
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
/**************************************************************/
/******************** Gap estimator Scheme ********************/
/**************************************************************/
    else if(Mode==3 || Mode==4){
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
            for(int i_sample=0;i_sample<Nsample;++i_sample){
                diagonal_update();
                construct_link_vertex_list();
                loop_update();
                flip_bit_operator();

                measure_with_propagate_state(i_sample);
            }
            estimator_fileout(Filename);

            if(it==2*Nit-1){
                if(Mode==3) gap_estimator_setup_workspace(L);
                else if(Mode==4){
                    int nk=3;
                    double wk[]={
                    PI,PI,PI,
                    PI*(1.0+2.0/4.0),PI,PI,
                    PI,PI*(1.0+2.0/8.0),PI};
                    momentum_setup_workspace(Nx,Ny,Nz,nk,wk,Nsample);
                }
                for(int k=0;k<Nblock;++k){
                    for(int i_sample=0;i_sample<Nsample;++i_sample){
                        diagonal_update();
                        construct_link_vertex_list();
                        loop_update();
                        flip_bit_operator();

                        measure_with_propagate_state(i_sample);
                        if(Mode==3){
                            gap_estimator_collect_data(Sequence,L,Sigma0,Sigmap,Nsite,Bond2index,StructFactor);
                            gap_estimator_calc_power_spectrum();
                        }
                        else if(Mode==4) momentum_collect_data(Sequence,L,Sigma0,Sigmap,Nsite,Bond2index,Beta,i_sample);
                    }
                    if(Mode==3) gap_estimator_ave_power_spectrum_fileout(Filename);
                    else if(Mode==4) momentum_calc_mean_fileout(Filename);
                    estimator_fileout(Filename);
                }
                if(Mode==3) gap_estimator_free_memory();
                else if(Mode==4) momentum_free_memory();
            }
            if(it%2==1){
                beta_doubling();
            }
        }
    }

    free_memory();
    return 0;
}
#endif
