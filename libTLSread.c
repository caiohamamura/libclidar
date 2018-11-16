#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLidVoxel.h"


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
/*write TLS point cloud from binary data*/

void writeTLSpointFromBin(char *namen,double *bounds,FILE *opoo)
{
  uint32_t i=0,nBeams=0;
  uint64_t buffSize=0;   /*buffer size*/
  uint64_t offset=0;
  double zen=0,az=0;
  double xCent=0,yCent=0,zCent=0;
  double x=0,y=0,z=0;
  uint32_t shotN=0;     /*shot number within this scan*/
  uint8_t nHits=0,j=0;  /*number of hits of this beam*/
  float r=0;            /*range*/
  float refl=0;         /*reflectance*/
  char *buffer=NULL;
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*skip to 4 bytes from the end*/
  if(fseek(ipoo,(long)-4,SEEK_END)){ 
    fprintf(stderr,"fseek error from end\n");
    exit(1);
  }
  if(fread(&nBeams,sizeof(uint32_t),1,ipoo)!=1){
    fprintf(stderr,"Error reading number of points\n");
    exit(1);
  }
  buffSize=ftell(ipoo);

  /*read data into buffer*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error to start\n");
    exit(1);
  }
  buffer=challoc((uint64_t)buffSize,"buffer",0);   /*allocate spave*/
  if(fread(&buffer[0],sizeof(char),buffSize,ipoo)!=buffSize){  /*read beams*/
    fprintf(stderr,"Error reading point data\n");
    exit(1);
  }
  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }


  /*loop along buffer and output*/
  offset=0;
  for(i=0;i<nBeams;i++){
    memcpy(&zen,&buffer[offset],8);
    zen*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&az,&buffer[offset],8);
    az*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&xCent,&buffer[offset],8);
    offset+=8;
    memcpy(&yCent,&buffer[offset],8);
    offset+=8;
    memcpy(&zCent,&buffer[offset],8);
    offset+=8;
    memcpy(&shotN,&buffer[offset],4);
    offset+=4;
    memcpy(&nHits,&buffer[offset],1);
    offset+=1;


    for(j=0;j<nHits;j++){
      memcpy(&r,&buffer[offset],4);
      offset+=4;
      memcpy(&refl,&buffer[offset],4);
      offset+=4;

      x=(double)r*sin(zen)*cos(az)+xCent;
      y=(double)r*sin(zen)*sin(az)+yCent;
      z=(double)r*cos(zen)+zCent;

      /*check bounds*/
      if((x>=bounds[0])&&(y>=bounds[1])&&(z>=bounds[2])&&(x<=bounds[3])&&(y<=bounds[4])&&(z<=bounds[5])){
        fprintf(opoo,"%.3f %.3f %.3f %f %d %d %f %f %f\n",x,y,z,refl,j,nHits,zen,az,r);
      }
    }/*hit loop*/

  }/*beam loop*/
  TIDY(buffer);
  return;
}/*writeTLSpointFromBin*/


/*##################################################################*/
/*read single TLS scan, all data*/

tlsScan *readTLSpolarBinary(char *namen)
{
  int i=0,offset=0;
  unsigned char j=0;
  uint64_t buffSize=0;   /*buffer size*/
  tlsScan *scan=NULL;
  double tempX=0,tempY=0,tempZ=0;
  char *buffer=NULL;
  FILE *ipoo=NULL;

  fprintf(stdout,"Reading %s ",namen);
  if(!(scan=(tlsScan *)calloc(1,sizeof(tlsScan)))){
    fprintf(stderr,"error scan allocation.\n");
    exit(1);
  }
  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }

  /*last 4 bytes state number of beams*/
  if(fseek(ipoo,(long)-4,SEEK_END)){ /*skip to 4 bytes from the end*/
    fprintf(stderr,"fseek error from end\n");
    exit(1);
  }
  if(fread(&scan->nBeams,sizeof(uint32_t),1,ipoo)!=1){
    fprintf(stderr,"Error reading number of points\n");
    exit(1);
  }
  fprintf(stdout,"There are %d TLS beams\n",scan->nBeams);

  /*determine file position to set size*/
  buffSize=ftell(ipoo);

  /*read data into buffer*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    fprintf(stderr,"fseek error to start\n");
    exit(1);
  }
  buffer=challoc((uint64_t)buffSize,"buffer",0);   /*allocate spave*/
  if(fread(&buffer[0],sizeof(char),buffSize,ipoo)!=buffSize){  /*read beams*/
    fprintf(stderr,"Error reading point data\n");
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*allocate structure space*/
  if(!(scan->beam=(tlsBeam *)calloc((long)scan->nBeams,sizeof(tlsBeam)))){
    fprintf(stderr,"error beam allocation. Allocating %d\n",scan->nBeams);
    exit(1);
  }

  /*load buffer into structure*/
  offset=0;
  scan->xOff=scan->yOff=scan->zOff=-1.0;
  for(i=0;i<scan->nBeams;i++){
    memcpy(&(scan->beam[i].zen),&buffer[offset],8);
    scan->beam[i].zen*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&(scan->beam[i].az),&buffer[offset],8);
    scan->beam[i].az*=M_PI/180.0; /*convert to radians*/
    offset+=8;
    memcpy(&tempX,&buffer[offset],8);
    offset+=8;
    memcpy(&tempY,&buffer[offset],8);
    offset+=8;
    memcpy(&tempZ,&buffer[offset],8);
    offset+=8;
    memcpy(&(scan->beam[i].shotN),&buffer[offset],4);
    offset+=4;
    memcpy(&(scan->beam[i].nHits),&buffer[offset],1);
    offset+=1;

    /*apply offset to coords*/
    if(scan->yOff<0.0){
      scan->xOff=tempX;
      scan->yOff=tempY;
      scan->zOff=tempZ;
    }

    scan->beam[i].x=(float)(tempX-scan->xOff);
    scan->beam[i].y=(float)(tempY-scan->yOff);
    scan->beam[i].z=(float)(tempZ-scan->zOff);

    if(scan->beam[i].nHits>0){
      scan->beam[i].r=falloc((int)scan->beam[i].nHits,"range",i);
      scan->beam[i].refl=falloc((int)scan->beam[i].nHits,"refl",i);
      for(j=0;j<scan->beam[i].nHits;j++){
        memcpy(&(scan->beam[i].r[j]),&buffer[offset],4);
        offset+=4;
        memcpy(&(scan->beam[i].refl[j]),&buffer[offset],4);
        offset+=4;
      }/*hit loop*/
    }/*hit check*/

  }/*load into array loop*/

  TIDY(buffer);
  return(scan);
}/*readTLSpolarBinary*/


/*##################################################################*/
/*read a single TLS scan within a voxel*/

tlsScan *readOneTLS(char *namen,voxStruct *vox,char useFracGap,tlsVoxMap *map,int fInd)
{
  int k=0,n=0;
  int pInd=0;
  int nBuff=0,vPlace=0;
  int xBin=0,yBin=0,zBin=0;
  int *voxList=NULL,nIntersect=0;
  int *markInt(int,int *,int);
  uint32_t *markUint32(int,uint32_t *,uint32_t);
  uint32_t j=0;
  float minR=0,maxR=0,lastHitR=0;
  double grad[3],*rangeList=NULL;
  double xCent=0,yCent=0,zCent=0;
  double x=0,y=0,z=0;
  tlsScan *scan=NULL,*tempTLS=NULL;
  char hasHit=0,doIt=0;

  /*max range of Riegl: OTHERS ARE SHORTER. COULD BE ADJUSTABLE*/
  maxR=300.0;

  /*allocate space*/
  if(!(scan=(tlsScan *)calloc(1,sizeof(tlsScan)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  if(map->mapFile==NULL){
    map->mapFile=iIalloc(vox->nVox,"voxel file map",0);        /*file per voxel*/
    if(!(map->mapPoint=(uint32_t **)calloc(vox->nVox,sizeof(uint32_t *)))){
      fprintf(stderr,"error in voxel point map allocation.\n");
      exit(1);
    }
    map->nIn=ialloc(vox->nVox,"voxel map number",0);        /*file per voxel*/
  }

  /*read all data into RAM*/
  tempTLS=readTLSpolarBinary(namen);

  nBuff=4*tempTLS->nBeams;
  if(!(scan->point=(tlsPoint *)calloc(nBuff,sizeof(tlsPoint)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  scan->nPoints=scan->nBeams=0;
  scan->beam=NULL;
  scan->xOff=tempTLS->xOff;
  scan->yOff=tempTLS->yOff;
  scan->zOff=tempTLS->zOff;
  /*fprintf(stdout,"Scan centre %f %f %f\n",scan->xOff,scan->yOff,scan->zOff);*/

  /*determine which are within bounds*/
  /*are we within 300 m of the bounds?*/
  if(((vox->bounds[0]-tempTLS->xOff)<=maxR)&&((vox->bounds[1]-tempTLS->yOff)<=maxR)&&\
     ((vox->bounds[2]-tempTLS->zOff)<=maxR)&&((vox->bounds[3]-tempTLS->xOff)>=(-1.0*maxR))&&\
     ((vox->bounds[4]-tempTLS->yOff)>=(-1.0*maxR))&&((vox->bounds[5]-tempTLS->zOff)>=(-1.0*maxR))){

    /*loop over beams*/
    for(j=0;j<tempTLS->nBeams;j++){
      /*avoid tilt mount if needed*/
      if(fabs(tempTLS->beam[j].zen)>=vox->maxZen)continue;  /*skip if zenith too high*/

      xCent=(double)tempTLS->beam[j].x+tempTLS->xOff;
      yCent=(double)tempTLS->beam[j].y+tempTLS->yOff;
      zCent=(double)tempTLS->beam[j].z+tempTLS->zOff;

      /*intersecting voxels*/
      grad[0]=tempTLS->beam[j].zen;
      grad[1]=tempTLS->beam[j].az;
      grad[2]=-99999.0;
      voxList=findVoxels(&(grad[0]),xCent,yCent,zCent,vox->bounds,\
                  &(vox->res[0]),&nIntersect,vox->nX,vox->nY,vox->nZ,&rangeList);

      if(nIntersect==0)continue;   /*if no voxels intersected*/

      /*gap fraction*/
      if(tempTLS->beam[j].nHits>0)lastHitR=tempTLS->beam[j].r[tempTLS->beam[j].nHits-1];
      else                        lastHitR=100000.0;

      for(k=0;k<nIntersect;k++){
        /*hits before voxel*/
        if(!useFracGap){  /*simple method. All hit until last return*/
          if(rangeList[k]<=lastHitR)vox->hits[fInd][voxList[k]]+=1.0;
          else                      vox->miss[fInd][voxList[k]]+=1.0;
        }else{            /*John's fractional method*/
          fprintf(stderr,"John's folly method not implemented yet\n");
          exit(1);
        }/*hits before voxel*/

        /*hits within voxel*/
        /*are we beyond the last return?*/
        doIt=1;
        if(k>0){
          minR=rangeList[k-1];
          if(rangeList[k-1]>lastHitR)doIt=0;  /*no information after this*/
        }else minR=0.0;

        /*add up total length of beams passing through*/
        vox->totVol[fInd][voxList[k]]+=rangeList[k]-minR;

        /*if not beyond, is it a hit or a miss in this voxel*/
        if(doIt){
          hasHit=0;
          for(n=0;n<tempTLS->beam[j].nHits;n++){
            if((tempTLS->beam[j].r[n]>=minR)&&(tempTLS->beam[j].r[n]<=rangeList[k])){
              hasHit=1;
              break;
            }
          }
          /*count up number of hits and misses within voxel*/
          if(hasHit)vox->inHit[fInd][voxList[k]]+=1.0;
          else      vox->inMiss[fInd][voxList[k]]+=1.0;
          /*count up volume sampled*/
          if(tempTLS->beam[j].r[n]<=lastHitR)vox->sampVol[fInd][voxList[k]]+=rangeList[k]-minR;  /*not last return*/
          else                               vox->sampVol[fInd][voxList[k]]+=tempTLS->beam[j].r[n]-minR;  /*last return*/
        }/*hits within voxel*/
      }/*voxel intersection loop*/
      /*record and map useful points*/
      if(vox->savePts){
        for(k=0;k<tempTLS->beam[j].nHits;k++){
          x=xCent+tempTLS->beam[j].r[k]*sin(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
          y=yCent+tempTLS->beam[j].r[k]*cos(tempTLS->beam[j].az)*sin(tempTLS->beam[j].zen);
          z=zCent+tempTLS->beam[j].r[k]*cos(tempTLS->beam[j].zen);

          /*check bounds and copy point if within*/
          if((x>=vox->bounds[0])&&(y>=vox->bounds[1])&&(z>=vox->bounds[2])&&\
             (x<=vox->bounds[3])&&(y<=vox->bounds[4])&&(z<=vox->bounds[5])){
            /*voxel space coordinates*/
            xBin=(int)((x-vox->bounds[0])/vox->res[0]);
            yBin=(int)((y-vox->bounds[1])/vox->res[1]);
            zBin=(int)((z-vox->bounds[2])/vox->res[2]);
            vPlace=zBin*vox->nX*vox->nY+yBin*vox->nX+xBin;

            /*are we within voxel space, to avoid rounding errors*/
            if((xBin<0)||(xBin>=vox->nX)||(yBin<0)||(yBin>=vox->nY)||(zBin<0)||(zBin>=vox->nZ)){
              continue;
            }

            /*mark TLS points*/
            scan->point[scan->nPoints].x=(float)(x-scan->xOff);  /*subtract offset to save disk space*/
            scan->point[scan->nPoints].y=(float)(y-scan->yOff);  /*subtract offset to save disk space*/
            scan->point[scan->nPoints].z=(float)(z-scan->zOff);  /*subtract offset to save disk space*/
            scan->point[scan->nPoints].r=tempTLS->beam[j].r[k];
            scan->point[scan->nPoints].refl=tempTLS->beam[j].refl[k];
            scan->point[scan->nPoints].hitN=k;
            scan->point[scan->nPoints].nHits=tempTLS->beam[j].nHits;
            /*map to voxels*/
fprintf(stdout,"nIn %d in %d\n",map->nIn[vPlace],vPlace);
            map->mapFile[vPlace]=markInt(map->nIn[vPlace],&(map->mapFile[vPlace][0]),fInd);
            map->mapPoint[vPlace]=markUint32(map->nIn[vPlace],&(map->mapPoint[vPlace][0]),scan->nPoints);
            map->nIn[vPlace]++;
            scan->nPoints++;
          }
        }/*hit loop*/
      }/*record point switch*/

      TIDY(rangeList);
      TIDY(voxList);
    }/*beam loop*/

    /*reallocate*/
    if((scan->nPoints>0)&&(scan->nPoints<tempTLS->nPoints)){
      if(!(scan->point=(tlsPoint *)realloc(scan->point,scan->nPoints*sizeof(tlsPoint)))){
        fprintf(stderr,"Error in reallocation, allocating %lu\n",scan->nPoints*sizeof(tlsPoint));
        exit(1);
      }
    }else if(scan->nPoints==0)TIDY(scan->point);

    /*determine gap fraction*/
    for(vPlace=0;vPlace<vox->nVox;vPlace++){
      for(k=0;k<map->nIn[vPlace];k++){
        pInd=map->mapPoint[vPlace][k];
        if((vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace])>0.0){
          scan->point[pInd].gap=vox->hits[fInd][vPlace]/(vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace]);
        }else scan->point[pInd].gap=1.0;
      }
    }/*gap fraction loop*/
  }/*voxel bound check*/

  /*tidy temporary space*/
  tempTLS=tidyTLScan(tempTLS);
  return(scan);
}/*readOneTLS*/


/*##################################################################*/
/*tidy multiple TLS structures*/

tlsScan *tidyTLScans(tlsScan *scans,int nScans)
{
  int i=0;

  if(scans){
    for(i=0;i<nScans;i++){
      if(scans[i].beam){
        for(i=0;i<scans[i].nBeams;i++){
          TIDY(scans[i].beam[i].refl);
          TIDY(scans[i].beam[i].r);
        }
        TIDY(scans[i].beam);
      }
      TIDY(scans[i].point);
    }
    TIDY(scans);
  }

  return(scans);
}/*tidyTLScans*/


/*##################################################################*/
/*tidy up TLS sctructure*/

tlsScan *tidyTLScan(tlsScan *scan)
{
  uint32_t i=0;

  if(scan){
    if(scan->beam){
      for(i=0;i<scan->nBeams;i++){
        TIDY(scan->beam[i].refl);
        TIDY(scan->beam[i].r);
      }
      TIDY(scan->beam);
    }
    TIDY(scan->point);
    TIDY(scan);
  }
  return(scan);
}/*tidyTLSscan*/


/*the end*/
/*##################################################################*/

