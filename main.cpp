#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "analyze.h"


int NX,NY,LEN,RUNS,JK_BIN_COUNT;
double KAPPA,EPSILON;
int STEPS,FRAMES;
//Jack Knife blocking variables
double topx[NXMAX][MAXRUN];
double topy[NXMAX][MAXRUN];
double bottomx[NXMAX][MAXRUN];
double bottomy[NXMAX][MAXRUN];

// Ribbon width 
double rwidth[NXMAX][MAXRUN];

int main(int argc, char **argv)
{
  FILE *pf,*vf,*rf;
  char init_strip[256],trajectory_file[256],profile_file[256],validrunfile[256],observable_file[1024],numvalidrun_file[256],runnum[256];
  int frame_cnt=0;
  double NORM;

   switch (argc){
     case 5:
       sscanf(argv[1],"%d",&NX);    
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&STEPS); 
       break;
     default:
       print_and_exit("Usage: %s NX NY KAPPA STEPS\n",
           argv[0]);
   }
 
  sprintf(validrunfile,"../Sim_dump_ribbon/L%d/W%d/k%.1f/valid_thermal_runs.log",NX,NY,KAPPA);
  sprintf(numvalidrun_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/numvalidruns.log",NX,NY,KAPPA);

  if(NULL==(vf=fopen(validrunfile,"r")))
	print_and_exit("I could not open file with simulation run numbers %s\n",validrunfile);

  if(NULL==(rf=fopen(numvalidrun_file,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",numvalidrun_file);

  fscanf(rf, "%s", runnum);
  RUNS = atoi(runnum);
  JK_BIN_COUNT=RUNS;
  //printf("NX = %d\n",NX);
  //printf("RUNS = %d\n",RUNS);

  FRAMES=STEPS/PERIOD;

  int initialize();

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",NX,NY);
  //printf("Init_strip.gsd : %s\n",init_strip);

  load_gsd(init_strip,0); 
  //printf("NX = %d\n",NX);
 
  //Jack Knife Blocking
  int iblo,r=0,frames,imes;
  for(iblo=0; iblo<NX; iblo++)//counter over the ribbon edge nodes
  {
	bottomx[iblo][RUNS]=0;
	bottomy[iblo][RUNS]=0;
	topx[iblo][RUNS]=0;
	topy[iblo][RUNS]=0;
	for(imes=0;imes<RUNS;imes++)
	{
		bottomx[iblo][imes]=0;
		bottomy[iblo][imes]=0;
		topx[iblo][imes]=0;
		topy[iblo][imes]=0;
	}
  }
	
  while (fscanf(vf, "%s", observable_file) == 1)// 1 is returned if fscanf reads a valid filepath
  {
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"%s/traj_thermal.gsd",observable_file);
	  //printf("Trajectory file being read : %s\n",trajectory_file);

	  //Looping through frames
	  for(frames=FRAMES/2;frames<FRAMES;frames++)
	  {
		load_gsd(trajectory_file,frames);
		//Data Blocking
		for(iblo=0; iblo<NX; iblo++)
                {
			bottomx[iblo][r]+=position[3*iblo];
			bottomy[iblo][r]+=position[3*iblo+1];
			topx[iblo][r]+=position[3*(N-NX+iblo)];
                        topy[iblo][r]+=position[3*(N-NX+iblo)+1];
                }
	  }
	  r++;//counter through the runs
	  
  }
  fclose(vf);

  for(iblo=0; iblo<NX; iblo++)//counter over the ribbon edge nodes
  {
	for(imes=0;imes<RUNS;imes++)
        {
		bottomx[iblo][imes]/=(FRAMES/2);
                bottomy[iblo][imes]/=(FRAMES/2);
                topx[iblo][imes]/=(FRAMES/2);
                topy[iblo][imes]/=(FRAMES/2);
        }
	for(imes=0;imes<RUNS;imes++)
        {
		bottomx[iblo][RUNS]+=bottomx[iblo][imes];
		bottomy[iblo][RUNS]+=bottomy[iblo][imes];
		topx[iblo][RUNS]+=topx[iblo][imes];
		topy[iblo][RUNS]+=topy[iblo][imes];
	}
	
	//Jack Knife Blocking
	for(imes=0;imes<RUNS;imes++)
        {
		bottomx[iblo][imes]=(bottomx[iblo][RUNS]-bottomx[iblo][imes])/(RUNS-1.);
		bottomy[iblo][imes]=(bottomy[iblo][RUNS]-bottomy[iblo][imes])/(RUNS-1.);
		topx[iblo][imes]=(topx[iblo][RUNS]-topx[iblo][imes])/(RUNS-1.);
		topy[iblo][imes]=(topy[iblo][RUNS]-topy[iblo][imes])/(RUNS-1.);
	}

	//Evaluate Ribbon width using Jack Knife blocked data
	for(imes=0;imes<RUNS;imes++)
        {
		rwidth[iblo][imes]=sqrt(pow((topx[iblo][imes]-bottomx[iblo][imes]),2)+pow((topy[iblo][imes]-bottomy[iblo][imes]),2));
  	}

	rwidth[iblo][RUNS]=0;
	for(imes=0;imes<RUNS;imes++)
        {
		rwidth[iblo][RUNS]+=rwidth[iblo][imes];
	}
	rwidth[iblo][RUNS]/=(RUNS*1.0);
   }
	
   //Normalizing using the left clamped width
   NORM = 1.0;//rwidth[0][RUNS];
   for(iblo=0; iblo<NX; iblo++)//counter over the ribbon edge nodes
   {
	rwidth[iblo][RUNS]/=NORM;
	//printf("%d\t%.8f\n",iblo,rwidth[iblo][RUNS]);
   }  

   //Jack knife error
       double jk_error[NXMAX],jk_error_term1[NXMAX],jk_error_term2[NXMAX];
    //          Jack Knife Error        
    for(int i=0;i<NX;i++)
    {
        jk_error[i]=0;
        jk_error_term1[i]=0;
        jk_error_term2[i]=0;
    }

    JK_BIN_COUNT=RUNS;
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term1[i] += rwidth[i][j] * rwidth[i][j];
        }
        jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
        
	for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term2[i] += (1.0/JK_BIN_COUNT) * rwidth[i][j];
        }
        jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        //      JK Error        
        jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
        printf ("%d\t%.8g\t%.8g\n",i,rwidth[i][RUNS],jk_error[i]);
    }


	  
  return 0;
}
