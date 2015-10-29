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

const int N=6 , T=50;
char seq_model[10000][51];
int num_seq_model[10000][50];
double alfa[50][6],beta[50][6],gama[50][6],xi[50][6][6];

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

static void load_seq_model(char seq[10000][51],int num_seq[10000][50] ,const char *filename)
{
    int i=0, j=0;
    FILE *fp;
    fp = fopen(filename, "r");
    for(i=0;i<10000;i++)
    {
       for(j=0;j<51;j++)
       {
        
            fscanf(fp,"%c",&seq[i][j]);
            switch(seq[i][j])
            {
                case 'A':
                    num_seq[i][j]=0;
                    break;
                case 'B':
                    num_seq[i][j]=1;
                    break;
                case 'C':
                    num_seq[i][j]=2;
                    break;
                case 'D':
                    num_seq[i][j]=3;
                    break;
                case 'E':
                    num_seq[i][j]=4;
                    break;    
                case 'F':
                    num_seq[i][j]=5;
                    break;
                default:
                    break;
            } 
                
        }    
   }


}

static void loadHMM( HMM *hmm, const char *filename )
{
   int i, j;
   FILE *fp = open_or_die( filename, "r");

   hmm->model_name = (char *)malloc( sizeof(char) * (strlen( filename)+1));
   strcpy( hmm->model_name, filename );

   char token[MAX_LINE] = "";
   while( fscanf( fp, "%s", token ) > 0 )
   {
      if( token[0] == '\0' || token[0] == '\n' ) continue;

      if( strcmp( token, "initial:" ) == 0 ){
         fscanf(fp, "%d", &hmm->state_num );

         for( i = 0 ; i < hmm->state_num ; i++ )
            fscanf(fp, "%lf", &( hmm->initial[i] ) );
      }
      else if( strcmp( token, "transition:" ) == 0 ){
         fscanf(fp, "%d", &hmm->state_num );

         for( i = 0 ; i < hmm->state_num ; i++ )
            for( j = 0 ; j < hmm->state_num ; j++ )
               fscanf(fp, "%lf", &( hmm->transition[i][j] ));
      }
      else if( strcmp( token, "observation:" ) == 0 ){
         fscanf(fp, "%d", &hmm->observ_num );

         for( i = 0 ; i < hmm->observ_num ; i++ )
            for( j = 0 ; j < hmm->state_num ; j++ )
               fscanf(fp, "%lf", &( hmm->observation[i][j]) );
      }
   }
}

static void dumpHMM( FILE *fp, HMM *hmm )
{
   int i, j;

   //fprintf( fp, "model name: %s\n", hmm->model_name );
   fprintf( fp, "initial: %d\n", hmm->state_num );
   for( i = 0 ; i < hmm->state_num - 1; i++ )
      fprintf( fp, "%.5lf ", hmm->initial[i]);
   fprintf(fp, "%.5lf\n", hmm->initial[ hmm->state_num - 1 ] );

   fprintf( fp, "\ntransition: %d\n", hmm->state_num );
   for( i = 0 ; i < hmm->state_num ; i++ ){
      for( j = 0 ; j < hmm->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", hmm->transition[i][j] );
      fprintf(fp,"%.5lf\n", hmm->transition[i][hmm->state_num - 1]);
   }

   fprintf( fp, "\nobservation: %d\n", hmm->observ_num );
   for( i = 0 ; i < hmm->observ_num ; i++ ){
      for( j = 0 ; j < hmm->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", hmm->observation[i][j] );
      fprintf(fp,"%.5lf\n", hmm->observation[i][hmm->state_num - 1]);
   }
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
static void cal_forward(HMM *hmm,int *seq)
{
    int i,j,k;
    for(i=0 ; i<hmm->observ_num;i++)
    {
       for(j=0;j<hmm->state_num;j++)
       {
           
            if(i==0)
                alfa[i][j]=hmm->initial[j] * hmm->observation[j][seq[i]];
            else
            {
               double tmp = 0;
               for(k=0 ; k<hmm->state_num;k++)
                   tmp += alfa[i-1][k] * hmm->transition[k][j];
                alfa[i][j]= tmp * hmm->observation[j][seq[i]];
            }
            
       }
    }
}

static void cal_backward(HMM *hmm,int *seq)
{
    int i,j,k;
    for(i=hmm->observ_num-1;i>=0;i--)
    {
        for(j=0;j<hmm->state_num;j++)
        {
            if(i==hmm->observ_num-1)
                beta[i][j]=1;
            else
            {
                double tmp = 0;
                for(k=0;k<hmm->state_num;k++)
                    tmp += hmm->transition[j][k] * hmm->observation[k][seq[i+1]] * beta[i+1][k];
                beta[i][j] = tmp ;
            }
        }
    }
}

static void cal_gama(HMM *hmm)
{
    int i,j;
    for(i=0;i<hmm->observ_num;i++)
    {
        double tmp =0;
        for(j=0;j<hmm->state_num;j++)
            tmp+=alfa[i][j]*beta[i][j];
        
        for(j=0;j<hmm->state_num;j++)
            gama[i][j] = alfa[i][j] * beta[i][j] / tmp;

        
    }
}

static void cal_xi(HMM *hmm,int *seq)
{
    int i,j,k;
    for(i=0;i<hmm->observ_num-1;i++) 
    {
        double tmp =0;
        for (j=0; j<hmm->state_num;j++)
        {
            for (k=0; k<hmm->state_num;k++)
            {
                tmp += alfa[i][j] * hmm->transition[j][k]* hmm->observation[k][seq[i+1]]*beta[i+1][k];
                 
            }
        }

        for (j=0; j<hmm->state_num; j++)
        {
            for (k=0; k<hmm->state_num; k++)
            {
                xi[i][j][k]=alfa[i][j]*hmm->transition[j][k]* hmm->observation[k][seq[i+1]]*beta[i+1][k] / tmp;     
            }
        }
    }
}

static void update_initial(HMM *hmm)
{
    int i;
    for(i=0;i<hmm->state_num;i++)
        hmm->initial[i] = gama[0][i];
}

static void update_transition(HMM *hmm)
{
    int i,j,k;
    for(i=0;i<hmm->state_num;i++)
    {
        double tmp1;
        for(j=0;j<hmm->observ_num-1;j++)
            tmp1 += gama[j][i];
        for(k=0;k<hmm->state_num;k++)
        {
            double tmp2;
            for(j=0;j<hmm->observ_num-1;j++)
                tmp2 += xi[j][i][k];
            hmm->transition[i][k]= tmp2/tmp1;
        }   
    }
}
static void update_observation(HMM *hmm,int *seq)
{
    int i,j,k;
    for(i=0;i<hmm->state_num;i++)
    {
        double tmp1[6]={0},tmp2=0;
        for(j=0;j<hmm->observ_num;j++)
        {
            tmp1[seq[j]] += gama[j][i];
            tmp2 += gama[j][i];
        }
        for(k=0;k<hmm->observ_num;k++)
            hmm->observation[i][k] = tmp1[k]/tmp2;
    }
}
#endif
