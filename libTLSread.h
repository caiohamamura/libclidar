
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
  uint32_t shotN; /*shot number within this scan*/
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
  uint32_t hitN;  /*hit number*/
  uint8_t nHits;  /*number of hits of this beam*/
}tlsPoint;


/*##########################################*/
/*TLS scan*/

typedef struct{
  tlsBeam *beam;     /*array of beams*/
  tlsPoint *point;   /*array of points*/
  double xOff;       /*offset to allow coords to be floats*/
  double yOff;       /*offset to allow coords to be floats*/
  double zOff;       /*offset to allow coords to be floats*/
  uint32_t nBeams;   /*number of beams in this scan*/
  uint32_t nPoints;  /*number of points in this scan*/
}tlsScan;


/*##########################################*/
/*function definitions*/

tlsScan *tidyTLScan(tlsScan *);
tlsScan *tidyTLScans(tlsScan *,int);
tlsScan *readTLSwithinVox(char **,int,voxStruct *,char,tlsVoxMap *);
tlsScan *readTLSpolarBinary(char *);


/*the end*/
/*##########################################*/

