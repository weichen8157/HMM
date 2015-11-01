#ifndef HMM_HEADER_
#define HMM_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef MAX_STATE
#	define MAX_STATE	10
#endif

#ifndef MAX_OBSERV
#	define MAX_OBSERV	26
#endif

#ifndef MAX_SEQ
#	define	MAX_SEQ		200
#endif

#ifndef MAX_LINE
#	define MAX_LINE 	256
#endif

typedef struct{
   char *model_name;
   int state_num;					//number of state
   int observ_num;					//number of observation
   double initial[MAX_STATE];			//initial prob.
   double transition[MAX_STATE][MAX_STATE];	//transition prob.
   double observation[MAX_OBSERV][MAX_STATE];	//observation prob.
} HMM;

static FILE *open_or_die( const char *filename, const char *ht )
{
   FILE *fp = fopen( filename, ht );
   if( fp == NULL ){
      perror( filename);
      exit(1);
   }

   return fp;
}


static int load_models( const char *listname, HMM *hmm, const int max_num )
{
   FILE *fp = open_or_die( listname, "r" );

   int count = 0;
   char filename[MAX_LINE] = "";
   while( fscanf(fp, "%s", filename) == 1 ){
      loadHMM( &hmm[count], filename );
      count ++;

      if( count >= max_num ){
         return count;
      }
   }
   fclose(fp);

   return count;
}

static void dump_models( HMM *hmm, const int num )
{
   int i = 0;
   for( ; i < num ; i++ ){ 
      //		FILE *fp = open_or_die( hmm[i].model_name, "w" );
      dumpHMM( stderr, &hmm[i] );
   }
}
