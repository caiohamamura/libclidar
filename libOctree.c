#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libOctree.h"

/*########################*/
/*# Functions for octrees #*/
/*# in lidar programs     #*/
/*#########################*/

/*#######################################*/
/*# Copyright 2006-2017, Steven Hancock #*/
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


/*#######################################*/
/*allocate space for top level*/

octTreeStruct *allocateOctree(int nLevels,int topN,double minX,double maxX,double minY,double maxY)
{
  int i=0;
  octTreeStruct *octree=NULL;
  float dX=0,dY=0;

  /*allocate space*/
  if(!(octree=(octTreeStruct *)calloc(1,sizeof(octTreeStruct)))){
    fprintf(stderr,"error octree allocation.\n");
    exit(1);
  }

  /*copy bounds*/
  octree->minX=minX;
  octree->maxX=maxX;
  octree->minY=minY;
  octree->maxY=maxY;

  /*determine resoplution*/
  dX=(float)(maxX-minX);
  dY=(float)(maxY-minY);
  octree->res=(dX>dY)?(int)(dX/(double)topN+1.0):(int)(dY/(double)topN+1.0);
  octree->nX=(int)(dX/octree->res);
  octree->nY=(int)(dY/octree->res);

  /*set pointers blank*/
  if(!(octree->tree=(treeStruct **)calloc(octree->nX*octree->nY,sizeof(treeStruct *)))){
    fprintf(stderr,"error octree allocation.\n");
    exit(1);
  }
  for(i=octree->nX*octree->nY-1;i>=0;i--)octree->tree[i]=NULL;
  octree->mapFile=NULL;
  octree->mapPoint=NULL;
  octree->nIn=NULL;

  return(octree);
}/*allocateOctree*/


/*#######################################*/
/*fill in octree*/

void fillOctree(double x,double y,double z,int nFile,uint32_t nPoint,octTreeStruct *octree)
{
  int xBin=0,yBin=0;
  double x0=0,y0=0;
  void mapOctree(int,octTreeStruct *,double,double,double,float,double,double);

  /*determine top level*/
  xBin=(int)((x-octree->minX)/(double)octree->res+0.5);
  yBin=(int)((y-octree->minY)/(double)octree->res+0.5);

  /*bounds check*/
  if((xBin>=0)&&(xBin<octree->nX)&&(yBin>=0)&&(yBin<octree->nY)){
    x0=(double)xBin*(double)octree->res+octree->minX;
    y0=(double)yBin*(double)octree->res+octree->minY;

    mapOctree(0,octree,x,y,z,octree->res,x0,y0);

  }/*bounds check*/

  return;
}/*fillOctree*/


/*#######################################*/
/*build octree map*/

void mapOctree(int level,octTreeStruct *octree,double x,double y,double z,float res,double x0,double y0)
{
  void mapOctree(int,octTreeStruct *,double,double,double,float,double,double);

  if(level<(octree->nLevel-1)){  /*keep recurssing*/
    mapOctree(level+1,octree,x,y,z,res/2.0,x0,y0);

  }else{  /*mark the points*/
    /*has this already been allocated*/

  }

  return;
}/*mapOctree*/


/*the end*/
/*#######################################*/

