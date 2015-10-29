#include "hmm.h"
#include <math.h>
#include <string.h>
int main(int argc, char *argv[])
{
/*
	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);
*/
    int i , j ;
	HMM hmm_initial;
	loadHMM( &hmm_initial, argv[1] );
	dumpHMM( stderr, &hmm_initial );
    load_seq_model(seq_model,num_seq_model,argv[2]);
    int len = strlen(argv[2]);
    char name[len-4];

    for(i=4;i<len;i++)
    {   
        name[i-4]=argv[2][i];
        //printf("%d",argv[2][i]);
    }
     
/*
    for(i=0;i<10000;i++)
    {
       for(j=0;j<50;j++)
       {
            printf("%d",num_seq_model[i][j]); 
       }   
    printf("\n"); 
    }
*/
    hmm_initial.state_num=N;                
    hmm_initial.observ_num=T; 
    //printf("state_num:%d\n", hmm_initial.state_num);
    //printf("observ_num:%d\n", hmm_initial.observ_num);
    for(i=0;i<10000;i++)
    {
        cal_forward(&hmm_initial,num_seq_model[i]);
        cal_backward(&hmm_initial,&num_seq_model[i]);
        cal_gama(&hmm_initial);
        cal_xi(&hmm_initial,&num_seq_model[i]);
        update_initial(&hmm_initial);
        update_transition(&hmm_initial);
        update_observation(&hmm_initial,&num_seq_model[i]); 
    }
/*
    for(i=0;i<50;i++)
    {
        for(j=0;j<6;j++)
        {
            printf("%d",alfa[i][j]);
        }
        printf("\n");
    }                                                             }
*/
	printf("%f\n", log(1.5) ); // make sure the math library is included
	return 0;
}
