#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libReadLVIS.h" 


/*#######################*/
/*# A library for       #*/
/*# handling LVIS files #*/
/*# S Hancock, 2017     #*/
/*#######################*/


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


/*########################################*/
/*check data sizes*/

void checkLVISsizes()
{
  if(sizeof(float)!=4){
    fprintf(stderr,"Size error\n");
    exit(1);
  }
  if(sizeof(double)!=8){
    fprintf(stderr,"Size error\n");
    exit(1);
  }
  if(sizeof(unsigned char)!=1){
    fprintf(stderr,"Size error\n");
    exit(1);
  }
  return;
}/*checkLVISsizes*/


/*########################################*/
/*read data*/

lvisData *readLVISdata(char *namen,int *nWaves)
{
  uint64_t i=0,len=0;
  uint64_t offset=0;
  lvisData *data=NULL;
  char *buffer=NULL;
  FILE *ipoo=NULL;


  /*open file*/
  if((ipoo=fopen(namen,"rb"))==NULL){
    fprintf(stderr,"Error opening input file %s\n",namen);
    exit(1);
  }


  /*read file size*/
  if(fseek(ipoo,(long)0,SEEK_END)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }
  len=ftell(ipoo);
  if(fseek(ipoo,(long)0,SEEK_SET)){
    fprintf(stderr,"fseek error\n");
    exit(1);
  }


  /*allocate space*/
  buffer=challoc(len,"buffer",0);
  *nWaves=(int)(len/sizeof(lvisData));
  if(!(data=(lvisData *)calloc(*nWaves,sizeof(lvisData)))){
    fprintf(stderr,"error data structure allocation.\n");
    exit(1);
  }

  /*read data*/
  if(fread(&(buffer[0]),sizeof(char),len,ipoo)!=len){
    fprintf(stderr,"error reading data\n");
    exit(1);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  /*copy data*/
  offset=0;
  for(i=0;i<(*nWaves);i++){
    memcpy(&(data[i].lfid),&(buffer[offset]),sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    memcpy(&(data[i].shotN),&(buffer[offset]),sizeof(uint32_t));
    offset+=sizeof(uint32_t);
    memcpy(&(data[i].az),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].zen),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].range),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].lvistime),&(buffer[offset]),sizeof(double));
    offset+=sizeof(double);
    memcpy(&(data[i].lon0),&(buffer[offset]),sizeof(double));
    offset+=sizeof(double);
    memcpy(&(data[i].lat0),&(buffer[offset]),sizeof(double));
    offset+=sizeof(double);
    memcpy(&(data[i].z0),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].lon431),&(buffer[offset]),sizeof(double));
    offset+=sizeof(double);
    memcpy(&(data[i].lat431),&(buffer[offset]),sizeof(double));
    offset+=sizeof(double);
    memcpy(&(data[i].z431),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].sigmean),&(buffer[offset]),sizeof(float));
    offset+=sizeof(float);
    memcpy(&(data[i].txwave[0]),&(buffer[offset]),sizeof(unsigned char)*80);
    offset+=sizeof(unsigned char)*80;
    memcpy(&(data[i].rxwave[0]),&(buffer[offset]),sizeof(unsigned char)*432);
    offset+=sizeof(unsigned char)*432;


    /*byteswap*/
    data[i].lfid=u32OneSwap(data[i].lfid);
    data[i].shotN=u32OneSwap(data[i].shotN);
    data[i].az=floOneSwap(data[i].az);
    data[i].zen=floOneSwap(data[i].zen);
    data[i].range=floOneSwap(data[i].range);
    data[i].lvistime=doOneSwap(data[i].lvistime);
    data[i].lon0=doOneSwap(data[i].lon0);
    data[i].lat0=doOneSwap(data[i].lat0);
    data[i].z0=floOneSwap(data[i].z0);
    data[i].lon431=doOneSwap(data[i].lon431);
    data[i].lat431=doOneSwap(data[i].lat431);
    data[i].z431=floOneSwap(data[i].z431);
    data[i].sigmean=floOneSwap(data[i].sigmean);
  }
  TIDY(buffer);

  return(data);
}/*readLVISdata*/


/*the end*/
/*####################################################*/

