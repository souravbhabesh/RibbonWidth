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
double topx[NXMAX][MAXRUN][MAXFRAMES];
double topy[NXMAX][MAXRUN][MAXFRAMES];
double bottomx[NXMAX][MAXRUN][MAXFRAMES];
double bottomy[NXMAX][MAXRUN][MAXFRAMES];

// Ribbon width 
double framewidth[NXMAX][MAXRUN][MAXFRAMES];
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
 
  //Initializing
  int iblo,r=0,frames,imes;
	
  while (fscanf(vf, "%s", observable_file) == 1)// 1 is returned if fscanf reads a valid filepath
  {
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"%s/traj_thermal.gsd",observable_file);
	  //printf("Trajectory file being read : %s\n",trajectory_file);

	  //Looping through frames
	  for(frames=FRAMES/2;frames<FRAMES;frames++)
	  {
		load_gsd(trajectory_file,frames);
		//Reading x and y coordinates of top and bottom boundary
		for(iblo=0; iblo<NX; iblo++)
                {
		  bottomx[iblo][r][frames-FRAMES/2]=position[3*iblo];
		  bottomy[iblo][r][frames-FRAMES/2]=position[3*iblo+1];
		  topx[iblo][r][frames-FRAMES/2]=position[3*(N-NX+iblo)];
		  topy[iblo][r][frames-FRAMES/2]=position[3*(N-NX+iblo)+1];
                }
	  }
	  r++;//counter through the runs
  }
  fclose(vf);

  //Evaluating the ribbon width at each frame in each run
  for(iblo=0;iblo<NX;iblo++)
  {
        for(imes=0;imes<RUNS;imes++)
        {
                for(frames=0;frames<FRAMES/2;frames++)
                {
                        framewidth[iblo][imes][frames]=sqrt(pow((topx[iblo][imes][frames]-bottomx[iblo][imes][frames]),2)+pow((topy[iblo][imes][frames]-bottomy[iblo][imes][frames]),2));
                }
        }
  }

  //Initializing width(each run) array
  for(iblo=0;iblo<NX;iblo++)
    {
        for(imes=0;imes<=RUNS;imes++)
        {
                rwidth[iblo][imes]=0;
        }
    }

  //Time Average of ribbon width in each run
  for(iblo=0;iblo<NX;iblo++)
  {
    for(imes=0;imes<RUNS;imes++)
    {
      for(frames=0;frames<FRAMES/2;frames++)
      {
        rwidth[iblo][imes]+=framewidth[iblo][imes][frames];
      }
      rwidth[iblo][imes]/=(FRAMES/2);//Time average
      rwidth[iblo][RUNS]+=rwidth[iblo][imes];
    }
  }

  //Jack Knife Blocking
  for(iblo=0;iblo<NX;iblo++)
  {
    for(imes=0;imes<RUNS;imes++)
    {
      rwidth[iblo][imes]=(rwidth[iblo][RUNS]-rwidth[iblo][imes])/(RUNS-1.);
    }
  }

  //Jack knife error
  double jk_error[NXMAX],jk_error_term1[NXMAX],jk_error_term2[NXMAX];

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
    //JK Error        
    jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
    printf ("%d\t%.8g\t%.8g\n",i,rwidth[i][RUNS]/(JK_BIN_COUNT*1.0),jk_error[i]);
 }
	  
  return 0;
}
