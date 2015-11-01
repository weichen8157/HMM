#include "hmm.h"
int main(int argc, char *argv[])
{
    if(argc!=4)
    {
        printf("Please input ITER model_int.txt seq_model model\n");
        exit(1);
    }
    int i , j ;

	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	dump_models( hmms, 5);

}
