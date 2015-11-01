#include "hmm.h"
int main(int argc, char *argv[])
{
/*
	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);
*/
    if(argc!=5)
    {
        printf("Please input ITER model_int.txt seq_model model\n");
        exit(1);
    }
    int i , j ;
printf("\noi\n");   
	HMM hmm_initial;
	loadHMM( &hmm_initial, argv[2] );
	dumpHMM( stderr, &hmm_initial );
printf("\ngo\n");
    load_seq_model(seq_model,num_seq_model,argv[3]);
printf("\nba\n");
/*
    for(i=0;i<50;i++)
    {
       for(j=0;j<6;j++)
       {
            printf("%d",alfa[i][j]); 
       }   
    printf("\n"); 
    }
*/
    
    //printf("state_num:%d\n", hmm_initial.state_num);
    //printf("observ_num:%d\n", hmm_initial.observ_num);
printf("\nya\n");    
    //for(i=0;i<argv[1];i++)
    //{
        for(j=0;j<4;j++)
        {
            cal_forward(&hmm_initial,num_seq_model[j]);
            cal_backward(&hmm_initial,num_seq_model[j]);
            cal_gama(&hmm_initial);
            cal_xi(&hmm_initial,num_seq_model[j]);
            update_initial(&hmm_initial);
            update_transition(&hmm_initial);
            update_observation(&hmm_initial,num_seq_model[j]);
        }
/*
        for(i=0;i<50;i++)
        {
            for(j=0;j<6;j++)
            {
             printf("%.5lf",gama[i][j]); 
            }   
        printf("\n"); 
        }
*/
    //}
printf("\nyaya\n");
    dumpHMM( stderr, &hmm_initial );
    output_model(&hmm_initial,argv[4]);
printf("\nyayaya\n");
	printf("%f\n", log(1.5) ); // make sure the math library is included
	return 0;
}
