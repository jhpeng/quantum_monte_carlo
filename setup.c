#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "scheme_diluted_bilayer_heisenberg.h"


int Help;
int Mode, Nx, Ny;
double P,Jbond,Beta,Beta_i,Beta_v,Beta_f;
int Nsample, Nther, Nblock, Nit;
int Seed;
char Filename[128];
void set_opt(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hx:y:j:b:n:k:t:s:f:m:p:i:v:f:z:"))!=-1){
        switch(c){
            case 'h':
                Help=1;
                printf("usage:\n");
                printf("\t-h print this help\n");
                printf("\t-m schem for the calculation\n");
                printf("\t\t 0 : diluted bilayer heisenberg model in finite temperature\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-j <J/Q ratio> default 1.0\n");
                printf("\t-b <beta> default 4.0\n");
                printf("\t-i <beta_i> default 1.0\n");
                printf("\t-f <beta_f> default 4.0\n");
                printf("\t-v <beta_v> default 0.2\n");
                printf("\t-n <Nsample> default 2000\n");
                printf("\t-k <Nblock> default 50\n");
                printf("\t-t <Nther> default 2000\n");
                printf("\t-e number of iteration for beta-doubing scheme\n");
                printf("\t-s <seed of random number generator> default 1\n");
                printf("\t-z <the file name of output data> default \"test.txt\"\n");
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
            case 'i':
                Beta_i=atof(optarg);
                break;
            case 'f':
                Beta_f=atof(optarg);
                break;
            case 'v':
                Beta_v=atof(optarg);
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
            case 'z':
                strcpy(Filename,optarg);
                break;
        }
    }
}

int main(int argc, char** argv){

    /*--------------default value----------------*/
    Beta = 4;
    Beta_i = 1;
    Beta_f = 4;
    Beta_v = 0.2;
    Seed = 9237912;
    Nx = 8;
    Ny = 8;
    Jbond = 1;
    P = 0.0;
    Nther = 2000;
    Nsample = 2000;
    Nblock = 50;
    Mode = 0;
    Nit = 5;

    /*----------------get option-----------------*/
    set_opt(argc,argv);
    if(Help) return 0;

    if(Mode==0) fine_temp_diluted_bilayer_heisenberg(Nx,Ny,P,Jbond,Beta_i,Beta_f,Beta_v,Nther,Nsample,Seed,Filename);
}
