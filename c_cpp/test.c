#include "hmm.h"
char test_data[2500][50];
int num_test_data[2500][50];
double delta[50][6];
int psi[50][6];
double model[5]={0};
int main(int argc, char *argv[])
{
    if(argc!=4)
    {
        printf("Please input modellist.txt testing_data.txt result.txt\n");
        exit(1);
    }
    int m,i,j,k,s ;
    
    HMM hmms[5];
    load_models( "modellist.txt", hmms, 5); 
    dump_models( hmms, 5);
    load_test_data(test_data,num_test_data,argv[2]);
    FILE* result=fopen(argv[3], "w");

for(s=0;s<2500;s++)
{
    for(m=0;m<5;m++)
    {
	for(i=0;i<T;i++)
	{
	    for(j=0;j<N;j++)
	    {
		if(i==0)
		    delta[i][j]=hmms[m].initial[j]*hmms[m].observation[num_test_data[s][i]][j];
		else
		{
		    double max=0;
		    for(k=0;k<N;k++)
		    {
			double z;
			z=delta[i-1][k]*hmms[m].transition[k][j];
			if(z>max)
			    max=z;
		    }
		    delta[i][j]=max*hmms[m].observation[num_test_data[s][i]][j];
		}
	    }
	}
	for(i=0;i<N;i++)
	{
	   if(delta[T-1][i]>model[m]) 
	       model[m]=delta[T-1][i];
	}
    }
    double max_model;
    int num_model=1;
    max_model=model[0];
    for(i=1;i<5;i++)
    {
	if(model[i]>max_model)
	{
	    max_model=model[i];
	    num_model=i+1;
	}
    }
    fprintf(result,"model_0%d.txt %g\n",num_model,max_model);
}
exit(0);
}
