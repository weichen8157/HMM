#include "hmm.h"
char test_data[2500][50];
int num_test_data[2500][50];
char real_answer[2500][13];
char final_answer[2500][13];
int psi[50][6];
double model[5];
int main(int argc, char *argv[])
{
    if(argc!=4)
    {
        printf("Please input modellist.txt testing_data.txt result.txt\n");
        exit(1);
    }
    int m,i,j,k,s;
    double accuracy=0,correct=0;
    HMM hmms[5];
    load_models( "modellist.txt", hmms, 5); 
    dump_models( hmms, 5);
    load_test_data(test_data,num_test_data,argv[2]);
    FILE* result=fopen(argv[3], "w");
    load_test_answer(real_answer,"testing_answer.txt");

    for(s=0;s<2500;s++)
    {
        model[0]=verterbi(&hmms[0],num_test_data[s]);
        model[1]=verterbi(&hmms[1],num_test_data[s]);
        model[2]=verterbi(&hmms[2],num_test_data[s]);
        model[3]=verterbi(&hmms[3],num_test_data[s]);
        model[4]=verterbi(&hmms[4],num_test_data[s]);
    
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

     	char tmp[10];
	itoa(num_model,tmp);
        strcpy(final_answer[s],"model_0"); 
	strcat(final_answer[s],tmp);
	strcat(final_answer[s],".txt\0");

	   
    }

    for(i=0;i<2500;i++)
    {
	int check=0;
	for(j=0;j<12;j++)
	{
	    if(final_answer[i][j]!=real_answer[i][j])
		check=1;
	}
	if(check==0)
	    correct++;
    }
    //if(strcmp(real_answer[i],final_answer[i])==0)
    /*
    for(i=0;i<2500;i++)
    {
	int check = 0;

        for(j=0;j<8;j++)
 	{
	    if(final_answer[i][j]!=real_answer[i][j])
		check=1;
	}
	if(check==0)
	    correct++;
    }
    */
    accuracy=correct/2500;   
    printf("correct:%.5lf accuracy:%.5lf\n",correct,accuracy);
    exit(0);
}
