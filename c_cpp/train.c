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
    double tmp1;
    unsigned int ITER;
    ITER = atoi(argv[1]);
    HMM hmm_initial;
    loadHMM( &hmm_initial, argv[2] );
    //dumpHMM( stderr, &hmm_initial );
    load_seq_model(seq_model,num_seq_model,argv[3]);
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
    
    for(i=0;i<ITER;i++)
    {
	int k,l;
	for(k=0;k<N;k++)
	{
		sum_gama_1[k]=0;
		sum_gama_2[k]=0;
		sum_gama_3[k]=0;
	}
	for(k=0;k<N;k++)
	{
	    for(l=0;l<N;l++)
   	    {
		sum_gama_o[k][l]=0;
		sum_xi[k][l]=0;
	    }
	}

        for(j=0;j<10000;j++)
        {
            cal_forward(&hmm_initial,num_seq_model[j]);
            cal_backward(&hmm_initial,num_seq_model[j]);
            cal_gama(&hmm_initial);
            cal_xi(&hmm_initial,num_seq_model[j]);
	    sum_of_gama_1();
	    sum_of_gama_2();
	    sum_of_gama_3(num_seq_model[j]);
	    sum_of_xi();
        }
	update_initial(&hmm_initial);
        update_transition(&hmm_initial);
        update_observation(&hmm_initial);

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
    }
    dumpHMM( stderr, &hmm_initial );
    output_model(&hmm_initial,argv[4]);
	//printf("%f\n", log(1.5) ); // make sure the math library is included
	return 0;
}
