#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "analyze.h"

double h_width[MAXFRAMES][NXMAX];
double h_bb[MAXFRAMES][NXMAX];


int initialize()
{
   for(int i=0;i<FRAMES;i++)
   {
        for(int j=0;j<2*NX;j++)
        {
           h_width[i][j]=0;
        }
        for(int k=0;k<NX;k++)
        {
           h_bb[i][k]=0;
        }
   }
   return 0;
}


/*      Backbone Height average         */
int bb_hgt(int frame)
{
 int j=0; //count over the backbone nodes
 h_bb[frame][0] = 0; //first two nodes are fixed 
 h_bb[frame][1] = 0;
 j = 2;
 for(int i=0;i<N;i++)
   {
        if(particleID[i]==4)
        {
                h_bb[frame][j] = position[3*i+2];
                j++;
        }
   }
 h_bb[frame][j] = 0;
 h_bb[frame][j+1] = 0;
/*
 if(frame == 0)
 {
        for(int i=0;i<NX;i++)
        {
                printf("%d\t%.8f\n",i,h_bb[frame][i]);
        }
 }
*/
 return 0;
}


/*      Width Height Average    */
int width_hgt(int frame)
{
  int k=0,j=0,k_cnt,j_cnt;
  for(int i=0;i<NX;i++)
  {
	do {
	   //if (frame==0)
		//printf("%d\t%d\t%d\t%.8f\n",N,i,(i/2)+2*k*NX,position[3*((i/2)+2*k*NX)+2]);         
	   h_width[frame][i] += position[3*(i+2*k*NX)+2];
	   k++;
	}while((i+2*k*NX) < N);
	k_cnt = k;
	k=0;
  }

  for(int i=0;i<NX;i++)
  {
        h_width[frame][i] = h_width[frame][i]/k_cnt;
        //if(frame == 0)
                //printf ("%d\t%.8f\t%d\t%d\n",i,h_width[frame][i],k_cnt,j_cnt);
  }
  return 0;
}


