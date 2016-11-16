#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libTLSread.h"


/*######################*/
/*# A library for      #*/
/*# handling tls files #*/
/*# S Hancock, 2016    #*/
/*######################*/


/*#######################################*/
/*# Copyright 2006-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of libCLidar.                                      #*/
/*#                                                                      #*/
/*# libCLidar is free software: you can redistribute it and/or modify    #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# libCLidar is distributed in the hope that it will be useful,         #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with libClidar.  If not, see <http://www.gnu.org/licenses/>.#*/
/*########################################################################*/


/*##################################################################*/
/*read TLS scans*/

tlsScan *readTLSscan(char *namen,char readBeam,char readPoint)
{
  tlsScan *scan=NULL;
  FILE *ipoo=NULL;

  /*initialise*/
  if(!(scan=(tlsScan *)calloc(1,sizeof(tlsScan)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  scan->beams=NULL;
  scan->points=NULL;
  scan->nBeams=0;
  scan->nPoints=0;


  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(scan);
}/*readTLSscan*/


/*##################################################################*/
/*tidy up TLS sctructure*/

tlsScan *tidyTLScan(tlsScan *scan)
{
  uint32_t i=0;

  if(scan){
    if(scan->beams){
      for(i=0;i<scan->nBeams;i++){
        TIDY(scan->beams[i].refl);
        TIDY(scan->beams[i].r);
      }
      TIDY(scan->beams);
    }
    TIDY(scan->points);
    TIDY(scan);
  }
  return(scan);
}/*tidyTLSscan*/


/*####################################*/
/*allocate tls voxel structure*/

tlsVoxStr *tlsVoxAllocate(int nFiles,float *vRes,double *bounds)
{
  int i=0,j=0;
  tlsVoxStr *vox=NULL;
  voxStruct *voxAllocate(int,float *,double *);


  if(!(vox=(tlsVoxStr *)calloc(1,sizeof(tlsVoxStr)))){
    fprintf(stderr,"error voxel structure allocation.\n");
    exit(1);
  }

  /*voxel structures*/
  vox->vox=voxAllocate(nFiles,&(vRes[0]),bounds);
  vox->toTLS=voxAllocate(nFiles,&(vRes[0]),bounds);

  /*mean factors for projecting area upwards*/
  vox->meanRefl=fFalloc(vox->vox->nScans,"meanRefl",0);
  vox->meanGap=fFalloc(vox->vox->nScans,"meanGap",0);
  vox->contN=iIalloc(vox->vox->nScans,"cal count",0);
  for(i=0;i<vox->vox->nScans;i++){
    vox->meanRefl[i]=falloc(vox->vox->nVox,"meanRefl",i+1);
    vox->meanGap[i]=falloc(vox->vox->nVox,"meanGap",i+1);
    vox->contN[i]=ialloc(vox->vox->nVox,"cal count",i+1);
    for(j=vox->vox->nVox-1;j>=0;j--){
      vox->meanRefl[i][j]=vox->meanGap[i][j]=0.0;
      vox->contN[i][j]=0;
    }
  }

  /*TLS specific things to map back to individual scans*/
  vox->mapFile=iIalloc(vox->vox->nVox,"voxel file map",0);
  vox->nIn=ialloc(vox->vox->nVox,"number in",0);
  if(!(vox->mapPoint=(uint32_t **)calloc(vox->vox->nVox,sizeof(uint32_t *)))){
    fprintf(stderr,"error voxel map allocation.\n");
    exit(1);
  }

  return(vox);
}/*tlsVoxAllocate*/

/*the end*/
/*##################################################################*/

