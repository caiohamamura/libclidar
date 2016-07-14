#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLidVoxel.h"


/*#########################*/
/*# Functions for dealing #*/
/*# with voxels           #*/
/*# S Hancock             #*/
/*# 6th November 2014     #*/
/*#########################*/

#define TOLERANCE 0.000001   /*tolerance for intersection tests*/
#define VTOL 0.01            /*tolerance for voxel finding*/

/*global to pass between functions*/
double tanZen=0,cosZen=0,sinZen=0;  /*to save calculations*/
double tanAz=0,cosAz=0,sinAz=0;


/*#######################################*/
/*voxels intersecting beam with width*/

int *beamVoxels(float *gradIn,double x0,double y0,double z0,double *bounds,double *res,int nX,int nY,int nZ,int *nPix,double beamRad)
{
  int i=0,j=0,k=0;
  int nAng=0;
  int tempPix=0;
  int *pixList=NULL;
  int *tempList=NULL;
  int *findVoxels(double *,double,double,double,double *,double *,int *,int,int,int,double **);
  int *markInt(int,int *,int);
  double grad[3];
  float ang=0,angStep=0;
  double x=0,y=0,z=0;
  char foundNew=0;

  nAng=90;
  angStep=2.0*M_PI/(float)nAng;

  /*central beam*/
  for(i=0;i<3;i++)grad[i]=(double)gradIn[i];
  pixList=findVoxels(&(grad[0]),x0,y0,z0,bounds,res,nPix,nX,nY,nZ,NULL);

  /*loop around rim of the beam*/
  for(i=0;i<nAng;i++){
    ang=(float)i*angStep;
    x=beamRad*sin(ang)+x0;  /*new start along edge of beam*/
    y=beamRad*cos(ang)+y0;  /*new start along edge of beam*/
    z=z0;     /*this should take into account the zenith angle of the beam*/

    /*find voxels intersected by the beam along that edge*/
    for(j=0;j<3;j++)grad[j]=(double)gradIn[j];
    tempList=findVoxels(&(grad[0]),x,y,z,bounds,res,&tempPix,nX,nY,nZ,NULL);
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
        (*nPix)++;
      } /*if new, mark it*/
    }/*temporary list loop*/
    TIDY(tempList);
  }/*sub step loop*/
  return(pixList);
}/*beamPixels*/


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
  double *coords=NULL,*iCoords=NULL,vThis[6];
  char angUp(double);
  char angRight(double,double);
  char angForward(double,double);
  char onEdge(double *,double *,double *,int);
  double *sideTest(double,double,double *,double *,double *);
  double *markDo(int,double *,double);
  void findClosestFacet(double *,double *,double,double);

  if(grad[2]>-9999.0){   /*grad is a Cartesian vector*/
    zen=atan2(sqrt(grad[0]*grad[0]+grad[1]*grad[1]),grad[2]);
    az=atan2(grad[0],grad[1]);
  }else{                 /*grad is a polar vector*/
    zen=grad[0];
    az=grad[1];
  }
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
  vCorn[1]=bounds[2];  /*minY*/
  vCorn[2]=bounds[4];  /*minZ*/
  vCorn[3]=bounds[1];  /*maxX*/
  vCorn[4]=bounds[3];  /*maxY*/
  vCorn[5]=bounds[5];  /*maxZ*/
  (*nPix)=0;


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
      iCoords=sideTest(zen,az,coords,vThis,vRes);
    }else if(angUp(zen)){  /*top plate intercept*/
      x=(vThis[5]-coords[2])*tanZen*sinAz+coords[0];
      y=(vThis[5]-coords[2])*tanZen*cosAz+coords[1];

      if((x>=vThis[0])&&(x<=vThis[3])&&(y>=vThis[1])&&(y<=vThis[4])){  /*through top*/
        iCoords=dalloc(3,"intersect coords",0);
        iCoords[0]=x;
        iCoords[1]=y;
        iCoords[2]=vThis[5];
      }else{   /*through a side*/
        iCoords=sideTest(zen,az,coords,vThis,vRes);
      }
    }else{                 /*bottom plate intercept*/
      x=(vThis[2]-coords[2])*tanZen*sinAz+coords[0];
      y=(vThis[2]-coords[2])*tanZen*cosAz+coords[1];

      if((x>=vThis[0])&&(x<=vThis[3])&&(y>=vThis[1])&&(y<=vThis[4])){  /*through bottom*/
        iCoords=dalloc(3,"intersect coords",0);
        iCoords[0]=x;
        iCoords[1]=y;
        iCoords[2]=vThis[2];
      }else{         /*through a side*/
        iCoords=sideTest(zen,az,coords,vThis,vRes);
      }
    }
    /*mark results*/
    if((xBin>=0)&&(xBin<vX)&&(yBin>=0)&&(yBin<vY)&&(zBin>=0)&&(zBin<vZ)){ /*bounds check*/
      pixList=markInt(*nPix,pixList,xBin+vX*yBin+vX*vY*zBin);
      if(rangeList)rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coords[0]-xCent)*\
        (coords[0]-xCent)+(coords[1]-yCent)*(coords[1]-yCent)+(coords[2]-zCent)*(coords[2]-zCent)));
      (*nPix)++;
    }/*bounds check*/

    /*fprintf(stdout,"Coords %g %g %g\n",coords[0],coords[1],coords[2]);*/
    TIDY(coords);    /*update coordinates*/
    coords=iCoords;
    iCoords=NULL;
  }/*voxel while loop*/

  /*mark the exit range too*/
  if(rangeList)rangeList[0]=markDo(*nPix,rangeList[0],sqrt((coords[0]-xCent)*\
     (coords[0]-xCent)+(coords[1]-yCent)*(coords[1]-yCent)+(coords[2]-zCent)*(coords[2]-zCent)));

  /*tidy arrays*/
  TIDY(iCoords);
  TIDY(coords);
  return(pixList);
}/*findVoxels*/


/*#######################################*/
/*test for side intersection*/

double *sideTest(double zen,double az,double *coords,double *vThis,double *vRes)
{
  double *iCoords=NULL;
  double d=0,xSep=0;
  double r=0;
  char angRight(double,double);

  iCoords=dalloc(3,"intersect coordinates",0);

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
  return(iCoords);
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
/*make a waveform from a point cloud image*/

void waveFromImage(char **rImage,float **wave,int numb,int rNx,int rNy)
{
  int i=0,j=0;
  char doneIt=0;

  for(i=0;i<numb;i++)wave[0][i]=wave[1][i]=0.0;


  /*turn range images into a waveforms*/
  for(i=rNx*rNy-1;i>=0;i--){ /*image x-y loop*/
    doneIt=0;
    for(j=0;j<numb;j++){ /*image bin loop*/
      if(rImage[j][i]>0){
        if(doneIt==0){
          wave[0][j]+=1.0; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/
          doneIt=1;
        }/*first hit only*/
        wave[1][j]+=1.0; /*appRefl*(rimRes*rimRes)/(M_PI*beamRad*beamRad);*/  /*hits per bin*/
      }
    }/*range image loop*/
  }/*range image bin loop*/
  return;
}/*waveFromImage*/

/*the end*/
/*###########################################################################*/

