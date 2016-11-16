
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


/*####################################*/
/*TLS beams, polar coords*/

typedef struct{
  double zen;     /*zenith*/
  double az;      /*azimuth*/
  float x;        /*beam origin*/
  float y;        /*beam origin*/
  float z;        /*beam origin*/
  uint32_t shotN; /*beam number*/  /*is this necessary? Could I save space without it*/
  uint8_t nHits;  /*number of hits of this beam*/
  float *r;       /*range*/
  float *refl;    /*reflectance*/
}tlsBeam;


/*##########################################*/
/*TLS point cloud*/

typedef struct{
  int bin;       /*bin number*/
  float x;       /*coordinate*/
  float y;       /*coordinate*/
  float z;       /*coordinate*/
  float gap;     /*voxel gap fraction*/
  float r;       /*range*/
  uint16_t refl; /*reflectance*/
}tlsPoint;


/*##########################################*/
/*TLS scan*/

typedef struct{
  tlsBeam *beams;   /*array of beams*/
  tlsPoint *points; /*array of points*/
  double xOff;       /*offset to allow coords to be floats*/
  double yOff;       /*offset to allow coords to be floats*/
  double zOff;       /*offset to allow coords to be floats*/
  uint32_t nBeams;   /*number of beams in this scan*/
  uint32_t nPoints;  /*number of points in this scan*/
}tlsScan;


/*####################################*/
/*voxel gap structure*/

typedef struct{
  /*the voxel*/
  voxStruct *vox;   /*TLS voxels. Gap within voxel*/
  voxStruct *toTLS; /*gap to TLS voxels*/
  float **meanGap;  /*mean minimum gap for voxels*/
  float **meanRefl; /*mean reflectance for voxels*/
  int **contN;

  /*TLS map*/
  int **mapFile;        /*file per voxel*/
  uint32_t **mapPoint;  /*point per voxel*/
  int *nIn;             /*number of points per voxel*/
}tlsVoxStr;


/*##########################################*/
/*function definitions*/

tlsScan *tidyTLScan(tlsScan *);
tlsScan *readTLSscan(char *,char,char);
tlsVoxStr *tlsVoxAllocate(int,float *,double *);
tlsVoxStr *tlsVoxAllocate(int,tlsVoxStr *)

/*the end*/
/*##########################################*/

