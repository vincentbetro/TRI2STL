#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "geometry.h"
#include "TRI2STL.h"
#include "Vector.h"
#include "journal.h"
#include "Spacing_Field.h"
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

//global variables
int num_procs = 1;
int my_rank = 0;
FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

int main(int argcs, char* pArgs[])
  {
  //declare variables (if not used as a library)
  int i, nf;
  const int bdim = 132; //buffer dim
  char **fnames; //geom file name storage
  char buff[bdim]; //buffer
  time_t tm;
  char *t_char;
  int nregions;
  int *regions = NULL;
    
  //create geometry (if not used as a library...normally already done)
  geometry* geom = new geometry();
  
  in_f = stdin;
  out_f = stdout;
  
  //create journal file
  if ((jou_f=fopen("TRI2STL.jou","w")) == NULL)
    {
    printf("\nCouldn't open TRI2STL journal file");
    exit(0);
    }
  
  //check for standard input (if not used as a library)
  if (--argcs < 1)
    {
      printf("\nNo input file specified!");
      printf("\nUsing standard input!");
    } 
  else
    {
    if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file <%s>\n",pArgs[argcs]);
        fprintf(stderr,"\n\nUsage: TRI2STL.XXX [batch_input_file]\n");
        exit(0);
      }
    }

  if (in_f != stdin)
    {
    sprintf(buff,"TRI2STL.out");

    if ((out_f=fopen(buff,"w")) == NULL)
     {
      fprintf(stderr,"\nCouldn't open file output file %s",buff);
      fclose(jou_f);
      exit(0);
      }
    }
  
  #ifdef _DEBUG
  buff[0] = '\0';
  
  sprintf(buff,"debug.%d",0);

  if ((debug_f=fopen(buff,"w")) == NULL)
  {
    fprintf(stderr,"\nCouldn't open debug output file %s",buff);
    fflush(stderr);
    fclose(jou_f);
    fclose(out_f);
    exit(0);
  }
  #endif

  //perform time check (if not used as a library)
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nTRI2STL run started at %s",t_char);

  //print program info 
  fprintf(out_f,"\n===========================================================================");
  fprintf(out_f,"\n COPYRIGHT 2003-2010 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA            ");
  fprintf(out_f,"\n                                                                           ");
  fprintf(out_f,"\n                    RIGHTS IN DATA                                         ");
  fprintf(out_f,"\n                                                                           ");
  fprintf(out_f,"\n THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT        ");
  fprintf(out_f,"\n    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO              ");
  fprintf(out_f,"\n         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL             ");
  fprintf(out_f,"\n              RIGHTS IN DATA CONTRACT CLAUSES.                             ");
  fprintf(out_f,"\n       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY          ");
  fprintf(out_f,"\n              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA                   ");
  fprintf(out_f,"\n                                                                           ");
  fprintf(out_f,"\n TRI2STL conversion code with TetGen built in                              ");
  fprintf(out_f,"\n NOTE: This data includes UT SimCenter at Chattanooga codes                ");
  fprintf(out_f,"\n that were developed under private non-government funding.                 ");
  fprintf(out_f,"\n This software is submitted with limited rights to use, reproduce,         ");
  fprintf(out_f,"\n and disclose this data for Government Purposes only.                      ");
  fprintf(out_f,"\n Requests for access to the software for non-governmental purposes         ");
  fprintf(out_f,"\n should be referrred to                                                    "); 
  fprintf(out_f,"\n                                                                           ");
  fprintf(out_f,"\n    Dr. Steve Karman                  Dr. Vincent Betro                    "); 
  fprintf(out_f,"\n    Steve-Karman@utc.edu              Vincent-Betro@utc.edu                "); 
  fprintf(out_f,"\n    423-425-5492  or  423-425-5470    423-425-5434                         "); 
  fprintf(out_f,"\n                                                                           ");
  fprintf(out_f,"\n    University of Tennessee SimCenter at Chattanooga                       "); 
  fprintf(out_f,"\n    701 East M. L. King Boulevard                                          "); 
  fprintf(out_f,"\n    Chattanooga, TN 37403                                                  "); 
  fprintf(out_f,"\n===========================================================================\n");
  fflush(out_f);

  //read number of geom files (if not used as a library...normally already done)
  journal(in_f, out_f, jou_f, "#Number of geometry files ->",nf);
  //check for adequate number of geom files
  if (nf <= 0)
    {
    fprintf(out_f,"\nNumber of geometry files must be > 0\n");
    fflush(out_f);
    exit(0);
    }

  //read geometry (if not used as a library...normally already done)
  fnames = (char**)malloc(nf*sizeof(char*));
  
  for (i=0; i < nf; i++)
    {
    fnames[i] = (char*)malloc(bdim*sizeof(char));
    sprintf(buff,"#File %d ->",i+1);
    journal(in_f, out_f, jou_f, buff, fnames[i]);
    }
    
  journal(in_f, out_f, jou_f,"#Enter number of regions to be saved from Tetgen ->", nregions);
  regions = new int[nregions];
  for (i = 0; i < nregions; i++)
    journal(in_f, out_f, jou_f,"#Enter number of region to be stitched ->", regions[i]);
    
  geom->read_geom(nf, fnames);
  
  //clean up fnames
  for (i=0; i < nf; i++)
    free(fnames[i]);
  free(fnames);
  
  //run p_opt as if it were a library with input parameters
  tri2stl_lib(geom, nregions, regions);
 
  //perform end time check
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nTRI2STL run completed at %s",t_char);
  fflush(out_f);

  delete geom;

  if (out_f != stdout)
    fclose(out_f);
    
  #ifdef _DEBUG
  if (in_f != stdin)
    fclose(debug_f);
  #endif
  
  return(0);
}

