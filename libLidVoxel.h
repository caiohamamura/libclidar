/*########################*/
/*# structures for voxel #*/
/*# lidar programs       #*/
/*########################*/

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



/*###################################################*/
/*voxels*/

typedef struct{
  int nScans;        /*number of contributing files*/
  float **hits;      /*hits per scan location*/
  float **miss;      /*misses per scan location*/
  float *rmse;       /*rmse of signal going in*/
  int *contN;        /*for normalising ALS*/
  int nVox;          /*total number of voxels*/
  int nX;            /*number of voxels*/
  int nY;            /*number of voxels*/
  int nZ;            /*number of voxels*/
  double res[3];     /*voxel resolution*/
  double bounds[6];  /*voxel bounds minX minY minZ maxX maxY maxZ*/
  char useRMSE;      /*switch to save RAM in voxelate*/
}voxStruct;


/*###################################################*/
/*lidar radiometric parameters for voxels*/

typedef struct{
  float minRefl;   /*minimum refletance value to scale between 0 and 1*/
  float maxRefl;   /*maximum refletance value to scale between 0 and 1*/
  float appRefl;   /*scale between TLS reflectance and size*/
  float beamTanDiv; /*tan of tls beam divergence*/
  float beamRad;    /*TLS start radius*/
  float minGap;    /*minimum gap fraction correction to apply*/
}lidVoxPar;


/*#####################################*/
/*structure to hold a range image*/

typedef struct{
  int nBins;        /*number of range bins*/
  int nX;           /*number of x bins*/
  int nY;           /*number of y bins*/
  float rRes;       /*range resolution*/
  float iRes;       /*image resolution*/
  double bounds[6]; /*minX minY minZ, maxX maxY maxZ*/
  double x0;        /*central coordinate*/
  double y0;        /*central coordinate*/
  double z0;        /*central coordinate*/
  char **image;     /*gap image at each range*/
  float grad[3];    /*vector along rage image*/
}rImageStruct;


/*####################################*/
/*voxel gap structure*/

typedef struct{
  /*the voxel*/
  voxStruct *vox;   /*TLS voxels. Gap within voxel*/
  voxStruct *toTLS; /*gap to TLS voxels*/
  float **meanGap;  /*mean minimum gap for voxels*/
  float **meanRefl; /*mean reflectance for voxels*/
  int **contN;      /*number of beams contributing to the above means*/

  /*TLS map*/
  int **mapFile;        /*file per voxel*/
  uint32_t **mapPoint;  /*point per voxel*/
  int *nIn;             /*number of points per voxel*/
}tlsVoxMap;

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
rImageStruct *allocateRangeImage(float,float,float,float *,double *,double *);
int *findVoxels(double *,double,double,double,double *,double *,int *,int,int,int,double **);
int *beamVoxels(float *,double,double,double,double *,double *,int,int,int,int *,double,double **,float);
voxStruct *voxAllocate(int,float *,double *,char);
voxStruct *tidyVox(voxStruct *);
void tidyVoxelMap(tlsVoxMap *,int);
void silhouetteImage(int,pCloudStruct **,tlsScan **,rImageStruct *,lidVoxPar *,int *,int,tlsVoxMap *);
void waveFromImage(rImageStruct *,float **,char,float);
void fillInRimageGround(rImageStruct *);
double *findVoxelBounds(int *,int,voxStruct *,tlsVoxMap *,float *,double,double,double);
void setWaveformRange(float *,double,float *,int,float);
void rotateX(double *,double);
void rotateZ(double *,double);
void readBoundsFromTLS(double *,char **,int);
void beamVoxelBounds(double *,float *,float,char,double *);

/*the end*/
/*###################################################*/

