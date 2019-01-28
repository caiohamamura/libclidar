#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLidVoxel.h"
#include "libLasProcess.h"


/*#########################*/
/*# Functions for dealing #*/
/*# with voxels           #*/
/*# S Hancock             #*/
/*# 6th November 2014     #*/
/*#########################*/


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



#define TOLERANCE 0.000001   /*tolerance for intersection tests*/
#define VTOL 0.01            /*tolerance for voxel finding*/

/*global to pass between functions*/
double tanZen=0,cosZen=0,sinZen=0;  /*to save calculations*/
double tanAz=0,cosAz=0,sinAz=0;


/*#############################################*/
/*make silhouette image from point cloud*/

void silhouetteImage(int nFiles,pCloudStruct **alsData,tlsScan *tlsData,rImageStruct *rImage,lidVoxPar *lidPar,int *voxList,int nIn,tlsVoxMap *map)
{
  int i=0,k=0,bin=0;
  int vInd=0,pInd=0,fInd=0;
  uint32_t j=0;
  float zen=0,az=0;
  double vect[3];
  void markPointSilhouette(double *,rImageStruct *,int,lidVoxPar *,float,uint16_t,double);

  /*angles for rotation*/
  zen=(float)atan2(sqrt((double)rImage->grad[0]*(double)rImage->grad[0]+(double)rImage->grad[1]*(double)rImage->grad[1]),(double)rImage->grad[2]);
  az=(float)atan2((double)rImage->grad[0],(double)rImage->grad[1]);


  if(alsData){   /*use ALS data*/
    for(i=0;i<nFiles;i++){
      for(j=0;j<alsData[i]->nPoints;j++){
        vect[0]=alsData[i]->x[j]-rImage->x0;
        vect[1]=alsData[i]->y[j]-rImage->y0;
        vect[2]=alsData[i]->z[j]-rImage->z0;
        /*rotate to x-y plane*/
        rotateZ(vect,(double)(-1.0*az));
        rotateX(vect,(double)(-1.0*zen));
        bin=(int)(vect[2]/rImage->rRes+0.5);

        if((bin>=0)&&(bin<rImage->nBins)){
          /*black out the points*/
          markPointSilhouette(&(vect[0]),rImage,bin,lidPar,alsData[i]->gap[j],alsData[i]->refl[j],0.0);
        }
      }/*point loop*/
    }/*file loop*/
  }else if(tlsData){   /*use TLS data*/
    for(k=0;k<nIn;k++){
      vInd=voxList[k];
      for(i=0;i<map->nIn[vInd];i++){
        fInd=map->mapFile[vInd][i];
        pInd=map->mapPoint[vInd][i];

        vect[0]=(double)tlsData[fInd].point[pInd].x+tlsData[fInd].xOff-rImage->x0;
        vect[1]=(double)tlsData[fInd].point[pInd].y+tlsData[fInd].yOff-rImage->y0;
        vect[2]=(double)tlsData[fInd].point[pInd].z+tlsData[fInd].zOff-rImage->z0;

        /*rotate to x-y plane*/
        rotateZ(vect,(double)(-1.0*az));
        rotateX(vect,(double)(-1.0*zen));
        bin=(int)(vect[2]/rImage->rRes+0.5);

        if((bin>=0)&&(bin<rImage->nBins)){
          /*black out the points*/
          markPointSilhouette(&(vect[0]),rImage,bin,lidPar,tlsData[fInd].point[pInd].gap,tlsData[fInd].point[pInd].refl,tlsData[fInd].point[pInd].r);
        } 
      }/*point in voxel loop*/
    }/*voxel loop*/
  }else{  /*no data. Something is wrong*/
    fprintf(stderr,"No data provided\n");
    exit(1);
  }

  return;
}/*silhouetteImage*/


/*############################################*/
/*mark lidar point in range image*/

void markPointSilhouette(double *coord,rImageStruct *rImage,int bin,lidVoxPar *lidPar,float gap,uint16_t refl,double r)
{
  int xInd=0,yInd=0,rPlace=0;
  int xStart=0,xEnd=0;
  int yStart=0,yEnd=0;
  int xIcent=0,yIcent=0;
  float rad=0;
  float maxRsepSq=0,rSepSq=0;

  if(gap<lidPar->minGap)gap=lidPar->minGap;
  rad=tlsPointSize(r,refl,lidPar->beamTanDiv,lidPar->beamRad,lidPar->minRefl,lidPar->maxRefl,lidPar->appRefl,gap); //)*lidPar->appRefl/gap;

  /*range image*/
  xIcent=(int)((coord[0]/(double)rImage->iRes)+0.5*(double)rImage->nX);
  yIcent=(int)((coord[1]/(double)rImage->iRes)+0.5*(double)rImage->nY);
  xStart=xIcent-(int)(rad/rImage->iRes+0.5);
  xEnd=xIcent+(int)(rad/rImage->iRes+0.5);
  yStart=yIcent-(int)(rad/rImage->iRes+0.5);
  yEnd=yIcent+(int)(rad/rImage->iRes+0.5);
  maxRsepSq=rad*rad;

  if(xStart<0)xStart=-1;      /*enforce bounds*/
  else if(xStart>=rImage->nX)xStart=rImage->nX;
  if(xEnd<0)xEnd=-1;      /*enforce bounds*/
  else if(xEnd>=rImage->nX)xEnd=rImage->nX;
  if(yStart<0)yStart=-1;
  else if(yStart>=rImage->nY)yStart=rImage->nY;
  if(yEnd<0)yEnd=-1;
  else if(yEnd>=rImage->nY)yEnd=rImage->nY;   /*enforce bounds*/

  /*mark centre point*/
  if((xIcent>=0)&&(xIcent<rImage->nX)&&(yIcent>=0)&&(yIcent<rImage->nY)){
    rPlace=yIcent*rImage->nX+xIcent;
    rImage->image[bin][rPlace]=1;
  }

  /*mark other points*/
  for(xInd=xStart;xInd<=xEnd;xInd++){
    if((xInd<0)||(xInd>=rImage->nX))continue;
    for(yInd=yStart;yInd<=yEnd;yInd++){
      if((yInd<0)||(yInd>=rImage->nY))continue;
      rSepSq=(float)((xInd-xIcent)*(xInd-xIcent)+(yInd-yIcent)*(yInd-yIcent))*rImage->iRes*rImage->iRes;
      if(rSepSq<=maxRsepSq){
        rPlace=yInd*rImage->nX+xInd;
        rImage->image[bin][rPlace]=1;
      }/*check within point*/
    }/*loop around point*/
  }/*loop around point*/

  return;
}/*markPointSilhouette*/


/*############################################*/
/*determine hit size*/

float tlsPointSize(double range,uint16_t refl,float tanDiv,float beamRad,float min,float max,float rhoApp,float gap)
{
  float d=0;
  float appRefl=0;
  float reflScale=0;

  appRefl=(float)max-(float)min;   /*scale from DN to size, takes phase func and albedo into account*/

  d=range*tanDiv+beamRad;                /*beam radius*/
  reflScale=((float)refl-(float)min)/appRefl;
  if(reflScale<0.0)     reflScale=0.0;   /*keep to bounds*/
  else if(reflScale>1.0)reflScale=1.0;   /*keep to bounds*/
  if(rhoApp<0.0)rhoApp=0.0;
  if(gap>TOLERANCE)d*=sqrt(reflScale*rhoApp/gap);   /*take optics into account*/
  else             d*=sqrt(reflScale*rhoApp/TOLERANCE);
  return(d);
}/*tlsPointSize*/


/*#############################################*/
/*allocate structure for range image*/

rImageStruct *allocateRangeImage(float beamRad,float rRes,float iRes,float *grad,double *origin,double *bounds)
{
  int i=0,k=0;
  float zen=0;
  rImageStruct *rImage=NULL;

  if(!(rImage=(rImageStruct *)calloc(1,sizeof(rImageStruct)))){
    fprintf(stderr,"error range image structure allocation.\n");
    exit(1);
  }

  rImage->x0=origin[0];
  rImage->y0=origin[1];
  rImage->z0=origin[2];
  if(grad){
    if(fabs(grad[0]+grad[1]+grad[2])>TOLERANCE){
      for(i=0;i<3;i++)rImage->grad[i]=grad[i];
    }else{
      rImage->grad[0]=rImage->grad[1]=0.0;
      rImage->grad[2]=-1.0;
    }
  }else{
    rImage->grad[0]=rImage->grad[1]=0.0;
    rImage->grad[2]=-1.0;
  }

  /*angles for rotation*/
  zen=(float)atan2(sqrt((double)rImage->grad[0]*(double)rImage->grad[0]+\
       (double)rImage->grad[1]*(double)rImage->grad[1]),(double)rImage->grad[2]);

  rImage->bounds[0]=-1.0*(double)beamRad;
  rImage->bounds[1]=-1.0*(double)beamRad;
  rImage->bounds[2]=0.0;
  rImage->bounds[3]=(double)beamRad;
  rImage->bounds[4]=(double)beamRad;
  rImage->bounds[5]=fabs(bounds[5]-bounds[2])*-1.0*(double)cos(zen);


  rImage->rRes=rRes;
  rImage->iRes=iRes;
  rImage->nBins=(int)((rImage->bounds[5]-rImage->bounds[2])/rImage->rRes+0.5);
  rImage->nX=(int)((rImage->bounds[3]-rImage->bounds[0])/rImage->iRes+0.5);
  rImage->nY=(int)((rImage->bounds[4]-rImage->bounds[1])/rImage->iRes+0.5);

  /*allocate image and set blank*/
  rImage->image=chChalloc(rImage->nBins,"range image",0);
  for(i=0;i<rImage->nBins;i++){
    rImage->image[i]=challoc((uint64_t)rImage->nX*(uint64_t)rImage->nY,"range image",i+1);
    for(k=rImage->nX*rImage->nY-1;k>=0;k--)rImage->image[i][k]=0;
  }
  return(rImage);
}/*allocateRangeImage*/


/*#############################################*/
/*add up hits and misses for a single beam*/

void countVoxGap(double x,double y,double z,float *grad,voxStruct *vox,int retNumb,int nRet,float beamRad,int numb)
{
  int i=0;
  int *voxList=NULL,nTot=0;
  double *rangeList=NULL;

  /*only do this for last returns per beam*/
  if(retNumb<nRet)return;

  /*check that a vector is there*/
  if(grad){
    if(fabs(grad[0]+grad[1]+grad[2])>TOLERANCE){
      /*determine which voxels are intersected*/
      voxList=beamVoxels(&(grad[0]),x,y,z,&(vox->bounds[0]),&(vox->res[0]),vox->nX,vox->nY,vox->nZ,&nTot,beamRad,&rangeList,-1.0);

      /*loop along intersected voxels*/
      for(i=0;i<nTot;i++){
        if(rangeList[i]<=0.0)vox->hits[numb][voxList[i]]+=1.0;
        else                 vox->miss[numb][voxList[i]]+=1.0;
      }/*intersecting voxel loop*/
      TIDY(rangeList);
      TIDY(voxList);
    }
  }

  return;
}/*countVoxGap*/


/*#######################################*/
/*voxels intersecting beam with width*/

int *beamVoxels(float *gradIn,double x0,double y0,double z0,double *bounds,double *res,int nX,int nY,int nZ,int *nPix,double beamRad,double **rangeList,float vRes)
{
  int i=0,j=0,k=0;
  int tempPix=0;
  int *pixList=NULL;
  int *tempList=NULL;
  int *findVoxels(double *,double,double,double,double *,double *,int *,int,int,int,double **);
  int *markInt(int,int *,int);
  double *markDo(int,double *,double);
  double *tempRange=NULL;
  double grad[3];
  float ang=0,angStep=0;  /*angular steps around edge of beam*/
  float rad=0,radRes=0;   /*radius to step along radial lines*/
  double x=0,y=0,z=0;
  char foundNew=0;

  /*determine angular resolution*/
  if(vRes>0.0)angStep=atan2(vRes/3.0,beamRad);
  else        angStep=2.0*M_PI/90.0;

  /*determine radial resolution*/
  if(vRes>=beamRad)radRes=beamRad;
  else             radRes=vRes/2.0;

  /*central beam*/
  for(i=0;i<3;i++)grad[i]=(double)gradIn[i];
  pixList=findVoxels(&(grad[0]),x0,y0,z0,bounds,res,nPix,nX,nY,nZ,rangeList);

  /*loop around rim of the beam*/
  ang=0.0;
  while(ang<2.0*M_PI){  /*angular loop*/
    rad=0.0;
    while(rad<=(float)beamRad){   /*radial loop*/
      x=rad*sin(ang)+x0;  /*new start along edge of beam*/
      y=rad*cos(ang)+y0;  /*new start along edge of beam*/
      z=z0;     /*this should take into account the zenith angle of the beam*/

      /*find voxels intersected by the beam along that edge*/
      for(j=0;j<3;j++)grad[j]=(double)gradIn[j];
      tempList=findVoxels(&(grad[0]),x,y,z,bounds,res,&tempPix,nX,nY,nZ,&tempRange);
      /*now sort through*/
      for(j=0;j<tempPix;j++){
        foundNew=1;
        for(k=0;k<(*nPix);k++){
          if(pixList[k]==tempList[j]){
            foundNew=0;
            break;
          }
        }/*final list loop*/
        if(foundNew==1){  /*if new, mark it*/
          pixList=markInt(*nPix,pixList,tempList[j]);
          rangeList[0]=markDo(*nPix,rangeList[0],tempRange[j]);
          (*nPix)++;
        } /*if new, mark it*/
      }/*temporary list loop*/
      TIDY(tempList);
      TIDY(tempRange);

      rad+=radRes;
    }/*radial loop*/
    ang+=angStep;
  }/*sub step loop*/

  return(pixList);
}/*beamVoxels*/


/*###########################################################################*/
/*find intersecting voxels*/

int *findVoxels(double *grad,double xCent,double yCent,double zCent,double *bounds,double *vRes,int *nPix,int vX,int vY,int vZ,double **rangeList)
{
  int *pixList=NULL;
  int *markInt(int,int *,int);
  int xBin=0,yBin=0,zBin=0;
  double zen=0,az=0;
  double vCorn[6];
  double x=0,y=0;
  double *coords=NULL,iCoords[3],vThis[6];
  char angUp(double);
  char angRight(double,double);
  char angForward(double,double);
  char onEdge(double *,double *,double *,int);
  double *markDo(int,double *,double);
  void sideTest(double,double,double *,double *,double *,double *);
  void findClosestFacet(double *,double *,double,double);

  if(grad[2]>-9999.0){   /*grad is a Cartesian vector*/
    zen=atan2(sqrt(grad[0]*grad[0]+grad[1]*grad[1]),grad[2]);
    az=atan2(grad[0],grad[1]);
  }else{                 /*grad is a polar vector*/
    zen=grad[0];
    az=grad[1];
  }

  (*nPix)=0;

  sinZen=sin(zen);
  cosZen=cos(zen);
  tanZen=tan(zen);
  sinAz=sin(az);
  cosAz=cos(az);
  tanAz=tan(az);

  coords=dalloc(3,"coords",0);
  coords[0]=xCent;
  coords[1]=yCent;
  coords[2]=zCent;

  vCorn[0]=bounds[0];  /*minX*/
  vCorn[1]=bounds[1];  /*minY*/
  vCorn[2]=bounds[2];  /*minZ*/
  vCorn[3]=bounds[3];  /*maxX*/
  vCorn[4]=bounds[4];  /*maxY*/
  vCorn[5]=bounds[5];  /*maxZ*/


  /*if we are outside test for intersection and move point*/
  if((coords[0]>vCorn[3])||(coords[1]>vCorn[4])||(coords[0]<vCorn[0])||\
     (coords[1]<vCorn[1])||(coords[2]>vCorn[5])||(coords[2]<vCorn[2])){
    /*for each voxel bound facet determine the range to. Reset coord as bound intersection*/
    findClosestFacet(&(coords[0]),&(vCorn[0]),zen,az);
  }/*outside but heading towards voxel space check*/


  while((coords[0]<=vCorn[3])&&(coords[1]<=vCorn[4])&&(coords[0]>=vCorn[0])&&\
        (coords[1]>=vCorn[1])&&(coords[2]<=vCorn[5])&&(coords[2]>=vCorn[2])){
    if((!angRight(zen,az))&&(onEdge(coords,vCorn,vRes,0)))  xBin=(int)((coords[0]-vCorn[0])/vRes[0])-1;
    else                                                    xBin=(int)((coords[0]-vCorn[0])/vRes[0]);
    if((!angForward(zen,az))&&(onEdge(coords,vCorn,vRes,1)))yBin=(int)((coords[1]-vCorn[1])/vRes[1])-1;
    else                                                    yBin=(int)((coords[1]-vCorn[1])/vRes[1]);
    if((!angUp(zen))&&(onEdge(coords,vCorn,vRes,2)))        zBin=(int)((coords[2]-vCorn[2])/vRes[2])-1;
    else                                                    zBin=(int)((coords[2]-vCorn[2])/vRes[2]);
    vThis[0]=(double)xBin*vRes[0]+vCorn[0];
    vThis[1]=(double)yBin*vRes[1]+vCorn[1];
    vThis[2]=(double)zBin*vRes[2]+vCorn[2];
    vThis[3]=vThis[0]+vRes[0];
    vThis[4]=vThis[1]+vRes[1];
    vThis[5]=vThis[2]+vRes[2];

    if((zen==(M_PI/2.0))||(zen==(-M_PI/2.0))){  /*a side*/
      sideTest(zen,az,coords,vThis,vRes,&(iCoords[0]));
    }else if(angUp(zen)){  /*top plate intercept*/
      x=(vThis[5]-coords[2])*tanZen*sinAz+coords[0];
      y=(vThis[5]-coords[2])*tanZen*cosAz+coords[1];

      if((x>=vThis[0])&&(x<=vThis[3])&&(y>=vThis[1])&&(y<=vThis[4])){  /*through top*/
        iCoords[0]=x;
        iCoords[1]=y;
        iCoords[2]=vThis[5];
      }else{   /*through a side*/
        sideTest(zen,az,coords,vThis,vRes,&(iCoords[0]));
      }
    }else{                 /*bottom plate intercept*/
      x=(vThis[2]-coords[2])*tanZen*sinAz+coords[0];
      y=(vThis[2]-coords[2])*tanZen*cosAz+coords[1];

      if((x>=vThis[0])&&(x<=vThis[3])&&(y>=vThis[1])&&(y<=vThis[4])){  /*through bottom*/
        iCoords[0]=x;
        iCoords[1]=y;
        iCoords[2]=vThis[2];
      }else{         /*through a side*/
        sideTest(zen,az,coords,vThis,vRes,&(iCoords[0]));
      }
    }

    /*catch rare cases of beams getting stuck between voxels. Needs resolving in a nicer way*/
    if((*nPix)>1000){
      fprintf(stdout,"nPix %d %f %f\n",*nPix,zen*180.0/M_PI,az*180.0/M_PI);
      (*nPix)=0;
      TIDY(coords);
      TIDY(pixList);
      TIDY(rangeList[0]);    
      return(pixList);
    }

    /*mark results*/
    if((xBin>=0)&&(xBin<vX)&&(yBin>=0)&&(yBin<vY)&&(zBin>=0)&&(zBin<vZ)){ /*bounds check*/
      pixList=markInt(*nPix,pixList,xBin+vX*yBin+vX*vY*zBin);
      if(rangeList)rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coords[0]-xCent)*\
        (coords[0]-xCent)+(coords[1]-yCent)*(coords[1]-yCent)+(coords[2]-zCent)*(coords[2]-zCent)));
      (*nPix)++;
    }/*bounds check*/

    coords[0]=iCoords[0];
    coords[1]=iCoords[1];
    coords[2]=iCoords[2];
  }/*voxel while loop*/

  /*mark the exit range too*/
  if(rangeList)rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coords[0]-xCent)*\
     (coords[0]-xCent)+(coords[1]-yCent)*(coords[1]-yCent)+(coords[2]-zCent)*(coords[2]-zCent)));

  /*tidy arrays*/
  TIDY(coords);
  return(pixList);
}/*findVoxels*/


/*#######################################*/
/*test for side intersection*/

void sideTest(double zen,double az,double *coords,double *vThis,double *vRes,double *iCoords)
{
  double d=0,xSep=0;
  double r=0;
  char angRight(double,double);

  if(angRight(zen,az))xSep=vRes[0]-(coords[0]-vThis[0]);
  else                xSep=vThis[0]-coords[0];
  d=xSep/tanAz;
  if((d+coords[1])>vThis[4]){        /*max Y*/
    iCoords[1]=vThis[4];
    r=(iCoords[1]-coords[1])/(sinZen*cosAz);
    iCoords[0]=r*sinZen*sinAz+coords[0];
    iCoords[2]=r*cosZen+coords[2];
  }else if((d+coords[1])<vThis[1]){  /*min Y*/
    iCoords[1]=vThis[1];
    r=(iCoords[1]-coords[1])/(sinZen*cosAz);
    iCoords[0]=r*sinZen*sinAz+coords[0];
    iCoords[2]=r*cosZen+coords[2];
  }else if(angRight(zen,az)){        /*max X*/
    iCoords[0]=vThis[3];
    r=(iCoords[0]-coords[0])/(sinZen*sinAz);
    iCoords[1]=r*sinZen*cosAz+coords[1];
    iCoords[2]=r*cosZen+coords[2];
  }else{                                /*min X*/
    iCoords[0]=vThis[0];
    r=(iCoords[0]-coords[0])/(sinZen*sinAz);
    iCoords[1]=r*sinZen*cosAz+coords[1];
    iCoords[2]=r*cosZen+coords[2];
  }
  return;
}/*sideTest*/


/*#######################################*/
/*are we on an edge*/

char onEdge(double *coords,double *vCorn,double *vRes,int dir)
{
  double tol=0;
  double test=0;

  tol=0.0001;
  test=(coords[dir]-vCorn[dir])/vRes[dir]-(double)(int)((coords[dir]-vCorn[dir])/vRes[dir]);
  if((test<tol)&&(test>-1.0*tol))return(1);
  else                           return(0);
}/*onEdge*/


/*#######################################*/
/*is it going right*/

char angRight(double zen,double az)
{
  if(zen<0.0)az+=M_PI;  /*mirror it*/

  if(((az>=0.0)&&(az<=M_PI))||((az>=2.0*M_PI)&&(az<=3.0*M_PI)))return(1);
  else                     return(0);
}/*angRight*/


/*#######################################*/
/*is it going forward*/

char angForward(double zen,double az)
{
  if(zen<0.0)az+=M_PI;  /*mirror it*/

  if(((az>=(-M_PI/2.0))&&(az<=(M_PI/2.0)))||((az>=(M_PI*3.0/2.0))&&(az<=(M_PI*5.0/2.0))))return(1);
  else                                                                                   return(0);
}/*angForward*/


/*#######################################*/
/*is it going up*/

char angUp(double zen)
{
  if(((zen>=(-M_PI/2.0))&&(zen<=(M_PI/2.0)))||((zen>=(M_PI*3.0/2.0))&&(zen<=(M_PI*5.0/2.0))))return(1);
  else                                                                                       return(0);
}/*angUp*/


/*#######################################*/
/*if outside bounds find closest facet*/

void findClosestFacet(double *coords,double *vCorn,double zen,double az)
{
  int i=0,j=0;
  double range[6];  /*range from centre to each facet*/
  double vect[3],thisC[3];
  double transCoord[3];
  double minR=0;
  char outside=0,found=0;

  /*lidar vector*/
  vect[0]=sin(zen)*sin(az);
  vect[1]=sin(zen)*cos(az);
  vect[2]=cos(zen);

  /*vCorn  minX minY minZ maxX maxY maxZ*/
  minR=100000.0;
  found=0;
  for(i=0;i<6;i++){/*facet loop*/
    range[i]=(vCorn[i]-coords[i%3])/vect[i%3];

    /*see if this is on the voxel*/
    outside=0;
    for(j=0;j<3;j++){
      thisC[j]=range[i]*vect[j]+coords[j];
      if((thisC[j]<vCorn[j])||(thisC[j]>vCorn[j+3])){
        outside=1;
        break;
      }
    }

    /*if not outside keep track of closest*/
    if(!outside){
      if((range[i]>=0.0)&&(range[i]<minR)){
        minR=range[i];
        for(j=0;j<3;j++)transCoord[j]=thisC[j];
        found=1;
      }
    }/*is on voxel facet check*/
  }/*facet loop*/

  /*copy if found*/
  if(found){
    for(i=0;i<3;i++)coords[i]=transCoord[i];
  }

  return;
}/*findClosestFacet*/


/*###############################################*/
/*make ground return slice solid*/

void fillInRimageGround(rImageStruct *rImage)
{
  int i=0,j=0,bin=0;
  char brEak=0;

  /*find lowest bin*/
  brEak=0;
  for(i=rImage->nBins-1;i>=0;i--){
    for(j=rImage->nX*rImage->nY-1;j>=0;j--){
      if(rImage->image[i][j]>0){
        brEak=1;
        bin=i;
        break;
      }
    }
    if(brEak)break;
  }
  if(brEak){
    for(j=rImage->nX*rImage->nY-1;j>=0;j--)rImage->image[bin][j]=1;
  }

  return;
}/*fillInRimageGround*/


/*###############################################*/
/*make a waveform from a point cloud image*/

void waveFromImage(rImageStruct *rImage,float **wave,char gaussFoot,float fSigma)
{
  int i=0,j=0,k=0;
  int place=0;
  float dx=0,dy=0;
  float weight=0,sep=0;
  float total=0,totWeight=0;
  char doneIt=0;

  for(i=0;i<rImage->nBins;i++){
    wave[0][i]=wave[1][i]=0.0;
  }
  if(gaussFoot==1){
    for(i=0;i<rImage->nX;i++){
      dx=(float)(i-rImage->nX/2)*rImage->iRes;
      for(j=0;j<rImage->nY;j++){
        dy=(float)(j-rImage->nY/2)*rImage->iRes;
        sep=sqrt(dx*dx+dy*dy);
        totWeight+=gaussian((double)sep,(double)fSigma,0.0);
      }
    }
  }else totWeight=(float)(rImage->nX*rImage->nY);

  /*turn range images into a waveforms*/
  for(i=0;i<rImage->nX;i++){
    dx=(float)(i-rImage->nX/2)*rImage->iRes;
    for(j=0;j<rImage->nY;j++){
      place=j*rImage->nX+i;
      dy=(float)(j-rImage->nY/2)*rImage->iRes;
      sep=sqrt(dx*dx+dy*dy);
      doneIt=0;

      for(k=0;k<rImage->nBins;k++){ /*image bin loop*/
        if(rImage->image[k][place]>0){
          if(gaussFoot==0){
            if(sep<=fSigma)weight=1.0;
            else           weight=0.0;
          }else if(gaussFoot==-1){
            weight=1.0;
          }else if(gaussFoot==1){
            weight=gaussian((double)sep,(double)fSigma,0.0);
          }

          if(doneIt==0){
            wave[0][k]+=weight; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/
            total+=weight;
            doneIt=1;
          }/*first hit only*/
          wave[1][k]+=weight; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/  /*hits per bin*/
        }
      }/*y loop*/
    }/*x loop*/
  }/*range image bin loop*/

  /*normalise waveforms*/
  total=0.0;
  for(j=0;j<rImage->nBins;j++)total+=wave[0][j];
  if(total>0.0){
    for(j=0;j<rImage->nBins;j++){
      wave[0][j]/=total;
      wave[1][j]/=totWeight;
    }
  }

  return;
}/*waveFromImage*/


/*############################################*/
/*allocate voxel structure*/

voxStruct *voxAllocate(int nFiles,float *vRes,double *bounds,char useRMSE)
{
  int i=0,j=0;
  voxStruct *vox=NULL;

  if(!(vox=(voxStruct *)calloc(1,sizeof(voxStruct)))){
    fprintf(stderr,"error voxel structure allocation.\n");
    exit(1);
  }

  /*note that findVoxels() needs minX maxX etc, different to dimage's minX minY etc*/
  for(i=0;i<3;i++)vox->res[i]=(double)vRes[i];
  for(i=0;i<6;i++)vox->bounds[i]=bounds[i];
  vox->nX=(int)((vox->bounds[3]-vox->bounds[0])/vox->res[0]+0.99);  /*add 0.99 to avoid rounding*/
  vox->nY=(int)((vox->bounds[4]-vox->bounds[1])/vox->res[1]+0.99);  /*add 0.99 to avoid rounding*/
  vox->nZ=(int)((vox->bounds[5]-vox->bounds[2])/vox->res[2]+0.99);  /*add 0.99 to avoid rounding*/
  vox->volume=(float)(vox->res[0]*vox->res[1]*vox->res[2]);

  vox->savePts=1;   /*defaults*/
  vox->maxZen=1000000.0;  /*use all points*/

  /*check for memory wrapping*/
  if(((uint64_t)vox->nX*(uint64_t)vox->nY*(uint64_t)vox->nZ)>=2147483647){
    fprintf(stderr,"Voxel bounds are too big to handle. Reduce %d %d %d\n",vox->nX,vox->nY,vox->nZ);
    exit(1);
  }

  vox->nVox=vox->nX*vox->nY*vox->nZ;
  vox->nScans=nFiles;

  vox->hits=fFalloc(vox->nScans,"voxel beam hits",0);
  vox->miss=fFalloc(vox->nScans,"voxel beam miss",0);
  vox->inHit=fFalloc(vox->nScans,"voxel point hits",0);
  vox->inMiss=fFalloc(vox->nScans,"voxel point miss",0);
  vox->sampVol=fFalloc(vox->nScans,"voxel volume sampled",0);
  vox->totVol=fFalloc(vox->nScans,"voxel volume total",0);
  vox->sumRsq=fFalloc(vox->nScans,"sum of radius of TLS points, squared",0);
  vox->meanRefl=fFalloc(vox->nScans,"mean reflectance of intersecting beams",0);
  vox->meanZen=fFalloc(vox->nScans,"mean zenith angle of intersecting beams",0);
  vox->contN=ialloc(vox->nVox,"voxel contribution",0);
  for(j=0;j<vox->nVox;j++)vox->contN[j]=0;

  vox->useRMSE=useRMSE;
  if(useRMSE){
    vox->rmse=falloc(vox->nVox,"voxel error",0);
    for(j=0;j<vox->nVox;j++)vox->rmse[j]=0.0;
  }

  for(i=0;i<vox->nScans;i++){
    vox->hits[i]=falloc(vox->nVox,"voxel hits",i+1);
    vox->miss[i]=falloc(vox->nVox,"voxel miss",i+1);
    vox->inHit[i]=falloc(vox->nVox,"voxel hits",i+1);
    vox->inMiss[i]=falloc(vox->nVox,"voxel miss",i+1);
    vox->sampVol[i]=falloc(vox->nVox,"voxel volume sampled",i+1);
    vox->totVol[i]=falloc(vox->nVox,"voxel volume total",i+1);
    vox->meanRefl[i]=falloc(vox->nVox,"mean reflectance of intersecting beams",i+1);
    vox->meanZen[i]=falloc(vox->nVox,"mean zenith angle of intersecting beams",i+1);
    vox->sumRsq[i]=falloc(vox->nVox,"sum of radius of TLS points, squared",i+1);
    for(j=0;j<vox->nVox;j++){
      vox->hits[i][j]=vox->miss[i][j]=vox->inHit[i][j]=vox->inMiss[i][j]=0.0;
      vox->sampVol[i][j]=vox->totVol[i][j]=vox->sumRsq[i][j]=0.0;
      vox->meanRefl[i][j]=vox->meanZen[i][j]=0.0;
    }
  }/*file loop*/

  return(vox);
}/*voxAllocate*/


/*#######################################*/
/*tidy voxel structure*/

voxStruct *tidyVox(voxStruct *vox)
{

  if(vox){
    TTIDY((void **)vox->hits,vox->nScans);
    TTIDY((void **)vox->miss,vox->nScans);
    TTIDY((void **)vox->inHit,vox->nScans);
    TTIDY((void **)vox->inMiss,vox->nScans);
    TTIDY((void **)vox->sampVol,vox->nScans);
    TTIDY((void **)vox->totVol,vox->nScans);
    TTIDY((void **)vox->meanRefl,vox->nScans);
    TTIDY((void **)vox->meanZen,vox->nScans);
    TTIDY((void **)vox->sumRsq,vox->nScans);
    TIDY(vox->rmse);
    TIDY(vox->contN);
    TIDY(vox);
  }

  return(NULL);
}/*tidyVox*/


/*###########################################################################*/
/*tidy up voxel map*/

void tidyVoxelMap(tlsVoxMap *map,int nVox)
{
  TTIDY((void **)map->mapFile,nVox);
  TTIDY((void **)map->mapPoint,nVox);
  TIDY(map->nIn);

  return;
}/*tidyVoxelMap*/


/*###########################################################################*/
/*find bounds of filled voxels*/

double *findVoxelBounds(int *voxList,int nIn,voxStruct *vox,tlsVoxMap *map,float *grad,double x0,double y0,double z0)
{
  int i=0,k=0,vInd=0;
  int ii=0,jj=0,kk=0;
  int xBin=0,yBin=0,zBin=0;
  float zen=0,az=0;
  double *bounds=NULL,vect[3];

  bounds=dalloc(6,"bounds",0);
  bounds[0]=bounds[1]=bounds[2]=10000000000.0;
  bounds[3]=bounds[4]=bounds[5]=-10000000000.0;

  zen=(float)atan2(sqrt((double)grad[0]*(double)grad[0]+(double)grad[1]*(double)grad[1]),(double)grad[2]);
  az=(float)atan2((double)grad[0],(double)grad[1]);

  /*loop over intersected voxels*/
  for(i=0;i<nIn;i++){
    vInd=voxList[i];
    if(map->nIn[vInd]>0){
      xBin=vInd%(vox->nX*vox->nY);
      yBin=(vInd-xBin)%vox->nX;
      zBin=((vInd-xBin)-yBin*vox->nX)/(vox->nX*vox->nY);

      /*each corner in turn*/
      for(ii=0;ii<2;ii++){
        for(jj=0;jj<2;jj++){
          for(kk=0;kk<2;kk++){
            /*rotate to beam vector*/
            vect[0]=(double)(xBin+ii)*vox->res[0]+vox->bounds[0]-x0;
            vect[1]=(double)(yBin+jj)*vox->res[1]+vox->bounds[1]-y0;
            vect[2]=(double)(zBin+kk)*vox->res[2]+vox->bounds[2]-z0;
            rotateZ(vect,(double)(-1.0*az));
            rotateX(vect,(double)(-1.0*zen));
            for(k=0;k<3;k++){   /*bound check*/
              if(vect[k]<bounds[k])bounds[k]=vect[k];
              if(vect[k]>bounds[k+3])bounds[k+3]=vect[k];
            }/*bound check*/
          }
        }
      }/*eight corner loop*/
    }/*filled voxel check*/
  }/*intersected voxel loop*/



  return(bounds);
}/*findVoxelBounds*/


/*###########################################################################*/
/*set waveform elevation along vector*/

void setWaveformRange(float *range,double z0,float *grad,int nBins,float res)
{
  int i=0;
  float r=0;
  double zen=0;

  zen=(float)atan2(sqrt((double)grad[0]*(double)grad[0]+(double)grad[1]*(double)grad[1]),(double)grad[2]);

  for(i=0;i<nBins;i++){
    r=(float)i*res;
    range[i]=z0+r*(float)cos(zen);
  }

  return;
}/*setWaveformRange*/


/*###########################################################################*/
/*read bounds for voxels from TLS*/

void readBoundsFromTLS(double *bounds,char **inList,int nScans)
{
  int i=0,k=0;
  uint32_t j=0,tInd=0;
  double x=0,y=0,z=0;
  double xCent=0,yCent=0,zCent=0;
  tlsScan *tempTLS=NULL;

  bounds[0]=bounds[1]=bounds[2]=10000000000.0;
  bounds[3]=bounds[4]=bounds[5]=-10000000000.0;

  for(i=0;i<nScans;i++){  /*file loop*/
    readTLSpolarBinary(inList[i],0,&tempTLS);
    for(j=0;j<tempTLS->nBeams;j++){/*point loop*/
      /*update TLS beams if needed*/
      readTLSpolarBinary(inList[i],j,&tempTLS);
      tInd=j-tempTLS->pOffset;   /*update index to account for buffered memory*/
      /*beam origin*/
      xCent=(double)tempTLS->beam[tInd].x+tempTLS->xOff;
      yCent=(double)tempTLS->beam[tInd].y+tempTLS->yOff;
      zCent=(double)tempTLS->beam[tInd].z+tempTLS->zOff;

      for(k=0;k<tempTLS->beam[tInd].nHits;k++){  /*hit loop*/
        /*point coordinate*/
        x=xCent+tempTLS->beam[tInd].r[k]*sin(tempTLS->beam[tInd].az)*sin(tempTLS->beam[tInd].zen);
        y=yCent+tempTLS->beam[tInd].r[k]*cos(tempTLS->beam[tInd].az)*sin(tempTLS->beam[tInd].zen);
        z=zCent+tempTLS->beam[tInd].r[k]*cos(tempTLS->beam[tInd].zen);

        /*determine bounds*/
        if(x<bounds[0])bounds[0]=x;
        if(y<bounds[1])bounds[1]=y;
        if(z<bounds[2])bounds[2]=z;
        if(x>bounds[3])bounds[3]=x;
        if(y>bounds[4])bounds[4]=y;
        if(z>bounds[5])bounds[5]=z;
      }/*hit loop*/
    }/*point loop*/
    tempTLS=tidyTLScan(tempTLS);
  }/*file loop*/

  return;
}/*readBoundsFromTLS*/


/*###########################################################################*/
/*clip x and y bounds to a beam*/

void beamVoxelBounds(double *origin,float *grad,float fSigma,char gaussFoot,double *bounds)
{
  int i=0;
  float rad=0;
  double x=0,y=0;

  if(gaussFoot)rad=determineGaussSep(fSigma,0.001);
  else         rad=fSigma;

  /*put beam start at top of voxel space*/
  origin[2]=bounds[5];

  bounds[0]=bounds[1]=100000000000.0;
  bounds[3]=bounds[4]=-100000000000.0;

  for(i=-1;i<=1;i+=2){  /*loop over edges*/
    /*top*/
    x=origin[0]+(float)i*rad;
    y=origin[1]+(float)i*rad;
    if(x<bounds[0])bounds[0]=x;
    if(x>bounds[3])bounds[3]=x;
    if(y<bounds[1])bounds[1]=y;
    if(y>bounds[4])bounds[4]=y;

    /*bottom*/
    x=origin[0]+(float)i*rad+((float)origin[2]-(float)bounds[2])*grad[0];
    y=origin[1]+(float)i*rad+((float)origin[2]-(float)bounds[2])*grad[1];
    if(x<bounds[0])bounds[0]=x;
    if(x>bounds[3])bounds[3]=x;
    if(y<bounds[1])bounds[1]=y;
    if(y>bounds[4])bounds[4]=y;
  }/*edges loop*/

  return;
}/*beamVoxelBounds*/


/*the end*/
/*###########################################################################*/

