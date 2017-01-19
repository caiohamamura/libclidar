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
  buffer=challoc(buffSize,"buffer",0);   /*allocate spave*/
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
/*read multiple TLS scan data within a voxel grid*/

tlsScan *readTLSwithinVox(char **inList,int nScans,voxStruct *vox,char useFracGap,tlsVoxMap *map)
{
  int i=0,k=0;
  int fInd=0,pInd=0;
  int nBuff=0,vPlace=0;
  int xBin=0,yBin=0,zBin=0;
  int *voxList=NULL,nIntersect=0;
  int *markInt(int,int *,int);
  uint32_t *markUint32(int,uint32_t *,uint32_t);
  uint32_t j=0;
  float maxR=0,lastHitR=0;
  double grad[3],*rangeList=NULL;
  double xCent=0,yCent=0,zCent=0;
  double x=0,y=0,z=0;
  tlsScan *scans=NULL,*tempTLS=NULL;

  /*max range of Riegl: OTHERS ARE SHORTER. COULD BE ADJUSTABLE*/
  maxR=300.0;

  /*allocate space*/
  if(!(scans=(tlsScan *)calloc(nScans,sizeof(tlsScan)))){
    fprintf(stderr,"error tls allocation.\n");
    exit(1);
  }
  map->mapFile=iIalloc(vox->nVox,"voxel file map",0);        /*file per voxel*/
  if(!(map->mapPoint=(uint32_t **)calloc(vox->nVox,sizeof(uint32_t *)))){
    fprintf(stderr,"error in voxel point map allocation.\n");
    exit(1);
  }
  map->nIn=ialloc(vox->nVox,"voxel map number",0);        /*file per voxel*/

  /*loop over scans*/
  for(i=0;i<nScans;i++){
    /*read all data into RAM*/
    tempTLS=readTLSpolarBinary(inList[i]);

    nBuff=4*tempTLS->nBeams;
    if(!(scans[i].point=(tlsPoint *)calloc(nBuff,sizeof(tlsPoint)))){
      fprintf(stderr,"error tls allocation.\n");
      exit(1);
    }
    scans[i].nPoints=scans[i].nBeams=0;
    scans[i].beam=NULL;
    scans[i].xOff=tempTLS->xOff;
    scans[i].yOff=tempTLS->yOff;
    scans[i].zOff=tempTLS->zOff;
    fprintf(stdout,"Scan centre %f %f %f\n",scans[i].xOff,scans[i].yOff,scans[i].zOff);

    /*determine which are within bounds*/
    /*are we within 300 m of the bounds?*/
    if(((vox->bounds[0]-tempTLS->xOff)<=maxR)&&((vox->bounds[1]-tempTLS->yOff)<=maxR)&&\
       ((vox->bounds[2]-tempTLS->zOff)<=maxR)&&((vox->bounds[3]-tempTLS->xOff)>=(-1.0*maxR))&&\
       ((vox->bounds[4]-tempTLS->yOff)>=(-1.0*maxR))&&((vox->bounds[5]-tempTLS->zOff)>=(-1.0*maxR))){


      /*loop over beams*/
      for(j=0;j<tempTLS->nBeams;j++){
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
          if(!useFracGap){  /*simple method. All hit until last return*/
            if(rangeList[k]<=lastHitR)vox->hits[i][voxList[k]]+=1.0;
            else                      vox->miss[i][voxList[k]]+=1.0;
          }else{            /*John's fractional method*/
            fprintf(stderr,"John's folly method not implemented yet\n");
            exit(1);
          }     
        }/*voxel intersection loop*/

        /*record and map useful points*/
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
            scans[i].point[scans[i].nPoints].x=(float)(x-scans[i].xOff);  /*subtract offset to save disk space*/
            scans[i].point[scans[i].nPoints].y=(float)(y-scans[i].yOff);  /*subtract offset to save disk space*/
            scans[i].point[scans[i].nPoints].z=(float)(z-scans[i].zOff);  /*subtract offset to save disk space*/
            scans[i].point[scans[i].nPoints].r=tempTLS->beam[j].r[k];
            scans[i].point[scans[i].nPoints].refl=tempTLS->beam[j].refl[k];
            scans[i].point[scans[i].nPoints].hitN=k;
            scans[i].point[scans[i].nPoints].nHits=tempTLS->beam[j].nHits;

            /*map to voxels*/
            map->mapFile[vPlace]=markInt(map->nIn[vPlace],&(map->mapFile[vPlace][0]),i);
            map->mapPoint[vPlace]=markUint32(map->nIn[vPlace],&(map->mapPoint[vPlace][0]),scans[i].nPoints);
            map->nIn[vPlace]++;

            scans[i].nPoints++;
          }
        }/*hit loop*/

        TIDY(rangeList);
        TIDY(voxList);
      }/*beam loop*/

      /*reallocate*/
      if((scans[i].nPoints>0)&&(scans[i].nPoints<tempTLS->nPoints)){
        if(!(scans[i].point=(tlsPoint *)realloc(scans[i].point,scans[i].nPoints*sizeof(tlsPoint)))){
          fprintf(stderr,"Balls\n");
          exit(1);
        }
      }else if(scans[i].nPoints==0)TIDY(scans[i].point);

      /*determine gap fraction*/
      for(vPlace=0;vPlace<vox->nVox;vPlace++){
        for(k=0;k<map->nIn[vPlace];k++){
          fInd=map->mapFile[vPlace][k];
          pInd=map->mapPoint[vPlace][k];
          if((vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace])>0.0){
            scans[fInd].point[pInd].gap=vox->hits[fInd][vPlace]/(vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace]);
          }else scans[fInd].point[pInd].gap=1.0;
        }
      }

    }/*voxel bound check*/

    tempTLS=tidyTLScan(tempTLS);
  }/*file loop*/

  return(scans);
}/*readTLSwithinVox*/


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

