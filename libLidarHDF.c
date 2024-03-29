#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "hdf5.h"
#include "stdint.h"
#include "tools.h"
#include "functionWrappers.h"
#include "libLidarHDF.h" 


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

int checkLVISsizes()
{
  if(sizeof(float)!=4){
    errorf("Size error\n");
    return(-1);
  }
  if(sizeof(double)!=8){
    errorf("Size error\n");
    return(-1);
  }
  if(sizeof(unsigned char)!=1){
    errorf("Size error\n");
    return(-1);
  }
  return(0);
}/*checkLVISsizes*/


/*########################################*/
/*read data*/

lvisLGWdata *readLVISlgw(char *namen,lvisLGWstruct *lvis)
{
  uint64_t i=0,len=0;
  uint64_t offset=0;
  uint16_t *swapUint16Arr(uint16_t *,int);
  float arg=0;
  lvisLGWdata *data=NULL;
  char *buffer=NULL;
  int lgwVersionFind(lvisLGWstruct *,char *,uint64_t);
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(namen,"rb"))==NULL){
    errorf("Error opening input file %s\n",namen);
    return(NULL);
  }

  /*read file size*/
  if(fseek(ipoo,(long)0,SEEK_END)){
    errorf("fseek error\n");
    return(NULL);
  }
  len=ftell(ipoo);
  if(fseek(ipoo,(long)0,SEEK_SET)){
    errorf("fseek error\n");
    return(NULL);
  }


  /*allocate reading space*/
  ASSIGN_CHECKNULL_RETNULL(buffer,challoc(len,"buffer",0));
  /*read data*/
  if(fread(&(buffer[0]),sizeof(char),len,ipoo)!=len){
    errorf("error reading data\n");
    return(NULL);
  }
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }

  /*determine version type*/
  ISINTRETNULL(lgwVersionFind(lvis,buffer,len));

  if(!(data=(lvisLGWdata *)calloc(lvis->nWaves,sizeof(lvisLGWdata)))){
    errorf("error data structure allocation.\n");
    return(NULL);
  }

  /*copy data*/
  offset=0;
  for(i=0;i<lvis->nWaves;i++){
    if(lvis->verMin>0){
      memcpy(&(data[i].lfid),&(buffer[offset]),sizeof(uint32_t));
      offset+=(uint64_t)sizeof(uint32_t);
      memcpy(&(data[i].shotN),&(buffer[offset]),sizeof(uint32_t));
      offset+=(uint64_t)sizeof(uint32_t);
    }
    if(lvis->verMin>=3){
      memcpy(&(data[i].az),&(buffer[offset]),sizeof(float));
      offset+=(uint64_t)sizeof(float);
      memcpy(&(data[i].zen),&(buffer[offset]),sizeof(float));
      offset+=(uint64_t)sizeof(float);
      memcpy(&(data[i].range),&(buffer[offset]),sizeof(float));
      offset+=(uint64_t)sizeof(float);
    }
    if(lvis->verMin>=2){
      memcpy(&(data[i].lvistime),&(buffer[offset]),sizeof(double));
      offset+=(uint64_t)sizeof(double);
    }
    memcpy(&(data[i].lon0),&(buffer[offset]),sizeof(double));
    offset+=(uint64_t)sizeof(double);
    memcpy(&(data[i].lat0),&(buffer[offset]),sizeof(double));
    offset+=(uint64_t)sizeof(double);
    memcpy(&(data[i].z0),&(buffer[offset]),sizeof(float));
    offset+=(uint64_t)sizeof(float);
    memcpy(&(data[i].lon431),&(buffer[offset]),sizeof(double));
    offset+=(uint64_t)sizeof(double);
    memcpy(&(data[i].lat431),&(buffer[offset]),sizeof(double));
    offset+=(uint64_t)sizeof(double);
    memcpy(&(data[i].z431),&(buffer[offset]),sizeof(float));
    offset+=(uint64_t)sizeof(float);
    memcpy(&(data[i].sigmean),&(buffer[offset]),sizeof(float));
    offset+=(uint64_t)sizeof(float);
    if(lvis->nTxBins>0){
      if(lvis->verMin<4){
        ASSIGN_CHECKNULL_RETNULL(data[i].txwave,uchalloc(lvis->nTxBins,"txwave",i+1));
        memcpy(&(data[i].txwave[0]),&(buffer[offset]),sizeof(unsigned char)*lvis->nTxBins);
        offset+=(uint64_t)sizeof(unsigned char)*(uint64_t)lvis->nTxBins;
      }else{   /*version 4 is a higher bit rate*/
        if(!(data[i].txwave4=(uint16_t *)calloc(lvis->nTxBins,sizeof(uint16_t)))){
          errorf("error in txwave allocation.\n");
          return(NULL);
        }
        memcpy(&(data[i].txwave4[0]),&(buffer[offset]),sizeof(uint16_t)*lvis->nTxBins);
        offset+=(uint64_t)sizeof(uint16_t)*(uint64_t)lvis->nTxBins;
      }
    }else data[i].txwave=NULL;
    if(lvis->verMin<4){
      ASSIGN_CHECKNULL_RETNULL(data[i].rxwave,uchalloc(lvis->nBins,"rxwave",0));
      memcpy(&(data[i].rxwave[0]),&(buffer[offset]),sizeof(unsigned char)*lvis->nBins);
      offset+=(uint64_t)sizeof(unsigned char)*(uint64_t)lvis->nBins;
    }else{   /*version 4 is a higher bit rate*/
      if(!(data[i].rxwave4=(uint16_t *)calloc(lvis->nBins,sizeof(uint16_t)))){
        errorf("error in rxwave allocation.\n");
        return(NULL);
      }
      memcpy(&(data[i].rxwave4[0]),&(buffer[offset]),sizeof(uint16_t)*lvis->nBins);
      offset+=(uint64_t)sizeof(uint16_t)*(uint64_t)lvis->nBins;
    }

    /*byteswap*/
    data[i].lfid=u32OneSwap(data[i].lfid);
    data[i].shotN=u32OneSwap(data[i].shotN);
    data[i].az=floOneSwap(data[i].az);
    data[i].range=floOneSwap(data[i].range);
    data[i].lvistime=doOneSwap(data[i].lvistime);
    data[i].lon0=doOneSwap(data[i].lon0);
    data[i].lat0=doOneSwap(data[i].lat0);
    data[i].z0=floOneSwap(data[i].z0);
    data[i].lon431=doOneSwap(data[i].lon431);
    data[i].lat431=doOneSwap(data[i].lat431);
    data[i].z431=floOneSwap(data[i].z431);
    data[i].sigmean=floOneSwap(data[i].sigmean);
    if(lvis->verMin>=3)data[i].zen=floOneSwap(data[i].zen);  /*only swap if read from file*/
    else{                                                    /*otherwise set from elevations*/
      arg=fabs(data[i].z0-data[i].z431)/(431.0*0.3);
      if(arg>0.0)data[i].zen=atan2(sqrt(1.0-arg*arg),arg)*180.0/M_PI;
    }
    if(lvis->verMin==4){   /*v1.4 needs byte swapping*/
      ASSIGN_CHECKNULL_RETNULL(data[i].rxwave4,swapUint16Arr(data[i].rxwave4,lvis->nBins));
    }
  }
  TIDY(buffer);

  return(data);
}/*readLVISlgw*/


/*#####################################*/
/*union for byte swapping*/

typedef union{
  char buff[sizeof(uint16_t)];
  uint16_t x;
}uint16Buff;  /*doubles*/


/*#####################################*/
/*byte swap a uint16_t array*/

uint16_t *swapUint16Arr(uint16_t *jimlad,int numb)
{
  int i=0,j=0;
  register int nBytes=sizeof(uint16_t);
  uint16_t *swap=NULL;
  uint16Buff ibuff,obuff;

  /*allocate space*/
  if(!(swap=(uint16_t *)calloc(numb,sizeof(uint16_t)))){
    errorf("error in uint16 swap array allocation.\n");
    return NULL;
  }

  /*loop over array*/
  for(i=0;i<numb;i++){
    ibuff.x=jimlad[i];
    /*sewap the bytes*/
    for(j=0;j<nBytes;j++)obuff.buff[j]=ibuff.buff[nBytes-1-j];
    swap[i]=obuff.x;
  }

  TIDY(jimlad);
  return(swap);
}/*swapUint16Arr*/


/*#####################################*/
/*determine lgw version*/

int lgwVersionFind(lvisLGWstruct *lvis,char *buffer,uint64_t len)
{
  int i=0,j=0,nWaves=0;
  int nVers=0,*rLen=NULL;
  uint64_t offset=0;
  double lat0=0;
  char thisVers=0;

  /*length of data packet in each version type*/
  nVers=5;
  ASSIGN_CHECKNULL_RETINT(rLen,ialloc(nVers,"record length",0));
  rLen[0]=476; /*(int)sizeof(struct lvis_lgw_v1_00);*/
  rLen[1]=484; /*(int)sizeof(struct lvis_lgw_v1_01);*/
  rLen[2]=492; /*(int)sizeof(struct lvis_lgw_v1_02);*/
  rLen[3]=584; /*(int)sizeof(struct lvis_lgw_v1_03);*/
  rLen[4]=1368; /*(int)sizeof(struct lvis_lgw_v1_04);*/
  lvis->verMin=0;
  lvis->verMaj=1;
  lvis->nWaves=0;

  /*do numbers make sense*/
  for(i=0;i<nVers;i++){
    nWaves=(int)(len/(uint64_t)rLen[i]);
    thisVers=1;
    for(j=0;j<nWaves;j++){
      offset=(uint64_t)j*(uint64_t)rLen[i];
      if(i==1)offset+=2*(uint64_t)sizeof(uint32_t)+(uint64_t)sizeof(double);
      else if(i==2)offset+=2*(uint64_t)sizeof(uint32_t)+2*(uint64_t)sizeof(double);
      else if(i>=3)offset+=2*(uint64_t)sizeof(uint32_t)+2*(uint64_t)sizeof(double)+3*(uint64_t)sizeof(float);

      memcpy(&lat0,&(buffer[offset]),sizeof(double));
      lat0=doOneSwap(lat0);

      if((lat0<-360.0)||(lat0>360.0)){
        thisVers=0;
        break;
      }//else msgf("poss %d %f %d\n",i,lat0,j);

    }/*wave loop*/

    /*is it this version?*/
    if(thisVers){
      lvis->verMin=i;
      lvis->nWaves=nWaves;
      msgf("LVIS version 1.%d\n",i);
      break;
    }
  }/*version loop*/

  /*set array lengths*/
  if(lvis->verMin<4)lvis->nBins=432;
  else              lvis->nBins=528;
  if(lvis->verMin<3)lvis->nTxBins=0;
  else if(lvis->verMin<4)lvis->nTxBins=80;
  else                   lvis->nTxBins=120;

  lvis->data=NULL;

  TIDY(rLen);
  return(0);
}/*lgwVersionFind*/


/*#####################################*/
/*tidy LVIS structure*/

lvisHDF *tidyLVISstruct(lvisHDF *lvis)
{
  if(lvis){
    TIDY(lvis->lon0);       /*LON0*/
    TIDY(lvis->lat0);       /*LAT0*/
    TIDY(lvis->lon1023);    /*LON1023*/
    TIDY(lvis->lat1023);    /*LAT1023*/
    TIDY(lvis->lfid);     /*LFID*/
    TIDY(lvis->shotN);    /*SHOTNUMBER*/
    TTIDY((void **)lvis->wave,1);    /*RXWAVE*/
    TTIDY((void **)lvis->pulse,1);   /*TXWAVE*/
    TIDY(lvis->zen);         /*INCIDENTANGLE*/
    TIDY(lvis->z0);         /*Z0*/
    TIDY(lvis->z1023);      /*Z1023*/
    TIDY(lvis->sigmean);     /*SIGMEAN*/
    TIDY(lvis->time);       /*TIME*/
    TIDY(lvis);
  }
  return(lvis);
}/*tidyLVISstruct*/


/*#####################################*/
/*readHDF5 LVIS file*/

lvisHDF *readLVIShdf(char *inNamen)
{
  int nWaves=0;
  lvisHDF *lvis=NULL;
  int checkNumber(int,int,char *);
  hid_t file;         /* Handles */
  char varName[10];  /*name for variable bin HDF5 files*/


  /*allocate structure*/
  if(!(lvis=(lvisHDF *)calloc(1,sizeof(lvisHDF)))){
    errorf("error in LVIS structure allocation.\n");
    return(NULL);
  }

  /*set to NULL to start with*/
  lvis->lon0=NULL;     /*LON0*/
  lvis->lat0=NULL;     /*LAT0*/
  lvis->lon1023=NULL;  /*LON1023*/
  lvis->lat1023=NULL;  /*LAT1023*/
  lvis->lfid=NULL;     /*LFID*/
  lvis->shotN=NULL;    /*SHOTNUMBER*/
  lvis->wave=NULL;     /*RXWAVE*/
  lvis->pulse=NULL;    /*TXWAVE*/
  lvis->zen=NULL;      /*INCIDENTANGLE*/
  lvis->z0=NULL;       /*Z0*/
  lvis->z1023=NULL;    /*Z1023*/
  lvis->sigmean=NULL;  /*SIGMEAN*/
  lvis->time=NULL;     /*TIME*/

  /*open HDF file*/
  msgf("Reading %s\n",inNamen);
  file=H5Fopen(inNamen,H5F_ACC_RDONLY,H5P_DEFAULT);

  /*read 1D double arrays*/
  ASSIGN_CHECKNULL_RETNULL(lvis->lon0,read1dDoubleHDF5(file,"LON0",&nWaves));
  lvis->nWaves=nWaves;
  ASSIGN_CHECKNULL_RETNULL(lvis->lat0,read1dDoubleHDF5(file,"LAT0",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"LAT0"));
  ASSIGN_CHECKNULL_RETNULL(lvis->lon1023,read1dDoubleHDF5(file,"LON1023",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"LON1023"));
  ASSIGN_CHECKNULL_RETNULL(lvis->lat1023,read1dDoubleHDF5(file,"LAT1023",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"LAT1023"));
  ASSIGN_CHECKNULL_RETNULL(lvis->time,read1dDoubleHDF5(file,"TIME",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"TIME"));

  /*read 1D float arrays*/
  ASSIGN_CHECKNULL_RETNULL(lvis->zen,read1dFloatHDF5(file,"INCIDENTANGLE",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"INCIDENTANGLE"));
  ASSIGN_CHECKNULL_RETNULL(lvis->z0,read1dFloatHDF5(file,"Z0",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"Z0"));
  ASSIGN_CHECKNULL_RETNULL(lvis->z1023,read1dFloatHDF5(file,"Z1023",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"Z1023"));
  ASSIGN_CHECKNULL_RETNULL(lvis->sigmean,read1dFloatHDF5(file,"SIGMEAN",&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"SIGMEAN"));

  /*read 1D uint32 arrays*/
  ASSIGN_CHECKNULL_RETNULL(lvis->lfid,read1dUint32HDF5(file,"LFID",&nWaves));
  NOT0RETNULL(checkNumber(nWaves,lvis->nWaves,"LFID"));
  ASSIGN_CHECKNULL_RETNULL(lvis->shotN,read1dUint32HDF5(file,"SHOTNUMBER",&nWaves));
  NOT0RETNULL(checkNumber(nWaves,lvis->nWaves,"SHOTNUMBER"));

  /*read 2d unit16 arrays*/
  ASSIGN_CHECKNULL_RETNULL(lvis->wave,read2dUint16HDF5(file,"RXWAVE",&lvis->nBins,&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,"RXWAVE"));
  /*if there is a pulse*/
  /*lvis->pulse=read2dUint16HDF5(file,"TXWAVE",&lvis->pBins,&nWaves);
  checkNumber(nWaves,lvis->nWaves,"TXWAVE");*/

  /*use the number of bins to determine the last bin variable names*/
  sprintf(varName,"LAT%d",lvis->nBins-1);
  ASSIGN_CHECKNULL_RETNULL(lvis->lat1023,read1dDoubleHDF5(file,varName,&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,varName));
  sprintf(varName,"LON%d",lvis->nBins-1);
  ASSIGN_CHECKNULL_RETNULL(lvis->lon1023,read1dDoubleHDF5(file,varName,&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,varName));
  sprintf(varName,"Z%d",lvis->nBins-1);
  ASSIGN_CHECKNULL_RETNULL(lvis->z1023,read1dFloatHDF5(file,varName,&nWaves));
  ISINTRETNULL(checkNumber(nWaves,lvis->nWaves,varName));

  /*close file*/
  if(H5Fclose(file)){
    errorf("Issue closing file\n");
    return(NULL);
  }

  return(lvis);
}/*readLVIShdf*/


/*#####################################*/
/*check integers match*/

int checkNumber(int newNumb,int oldNumb,char *label)
{
  if(newNumb!=oldNumb){
    errorf("Number mismatch %d %d for %s\n",newNumb,oldNumb,label);
    return(-1);
  }
  return(0);
}/*checkNumber*/


/*#####################################*/
/*read a HDF5 dataset*/

uint16_t **read2dUint16HDF5(hid_t file,char *label,int *nBins,int *nWaves)
{
  int i=0,ndims=0;
  uint16_t **jimlad=NULL;
  hid_t dset,space;
  hsize_t *dims=NULL;

  /*open dataset*/
  dset=H5Dopen2(file,label,H5P_DEFAULT);

  /*get dimensions*/
  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_ndims(space);
  if(!(dims=(hsize_t *)calloc(ndims,sizeof(hsize_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }

  if(H5Sget_simple_extent_dims(space,dims,NULL)!=ndims){
    errorf("Error\n");
    return(NULL);
  }
  (*nWaves)=(int)dims[0];
  (*nBins)=(int)dims[1];

  //tid=H5Dget_type(dset);

  /*allocate space*/
  jimlad=(uint16_t **)malloc(dims[0]*sizeof(uint16_t *));
  jimlad[0]=(uint16_t *)malloc(dims[0]*dims[1]*sizeof(uint16_t));
  for(i=1;i<dims[0];i++)jimlad[i]=jimlad[0]+i*dims[1];

  /*read data*/
  if(H5Dread(dset,H5T_NATIVE_USHORT,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad[0])){
    errorf("Error reading data %s\n",label);
    return(NULL);
  }

  /*close dataset*/
  if(H5Dclose(dset)){
    errorf("Error closing data %s\n",label);
    return(NULL);
  }
  if(H5Sclose(space)){
    errorf("Error closing space %s\n",label);
    return(NULL);
  }

  TIDY(dims);
  return(jimlad);
}/*read2dUint16HDF5*/


/*#########################################################*/
/*read a 1.5D float HDF5 dataset (2D compressed into 1D)*/

float *read15dFloatHDF5(hid_t file,char *label,int *nWaves,int *nBins)
{
  int ndims=0;
  float *jimlad=NULL;
  hid_t dset,space,filetype;
  hsize_t *dims=NULL;

  /*open dataset*/
  dset=H5Dopen2(file,label,H5P_DEFAULT);

  /*get dimensions*/
  space=H5Dget_space(dset);
  filetype=H5Dget_type(dset);
  ndims=H5Sget_simple_extent_ndims(space);
  if(!(dims=(hsize_t *)calloc(ndims,sizeof(hsize_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }

  if(H5Sget_simple_extent_dims(space,dims,NULL)!=ndims){
    errorf("Error\n");
    return(NULL);
  }
  (*nWaves)=(int)dims[0];
  (*nBins)=(int)dims[1];

  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(jimlad,falloc((uint64_t)(*nWaves)*(uint64_t)(*nBins),"1.5d float array HDF",0));

  /*read data*/
  if(H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad)){
    errorf("Error reading data %s\n",label);
    return(NULL);
  }

  /*close dataset*/
  if(H5Dclose(dset)){
    errorf("Error closing data %s\n",label);
    return(NULL);
  }
  if(H5Sclose(space)){
    errorf("Error closing space %s\n",label);
    return(NULL);
  }
  if(H5Tclose(filetype)){
    errorf("Error closing type %s\n",label);
    return(NULL);
  }

  TIDY(dims);
  return(jimlad);
}/*read15dFloatHDF5*/



/*#########################################################*/
/*read a 1.5D char HDF5 dataset (2D compressed into 1D)*/

char *read15dCharHDF5(hid_t file,char *label,int *nWaves,int *nBins)
{ 
  int ndims=0;
  char *jimlad=NULL;
  hid_t dset,space,filetype;
  hsize_t *dims=NULL;
  
  /*open dataset*/
  dset=H5Dopen2(file,label,H5P_DEFAULT);
  
  /*get dimensions*/
  space=H5Dget_space(dset);
  filetype=H5Dget_type(dset);
  ndims=H5Sget_simple_extent_ndims(space);
  if(!(dims=(hsize_t *)calloc(ndims,sizeof(hsize_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }
  
  if(H5Sget_simple_extent_dims(space,dims,NULL)!=ndims){
    errorf("Error\n");
    return(NULL);
  }
  (*nWaves)=(int)dims[0];
  (*nBins)=(int)dims[1];
  
  /*allocate space*/
  ASSIGN_CHECKNULL_RETNULL(jimlad,challoc((*nWaves)*(*nBins),"1.5d float array HDF",0));
  
  /*read data*/
  if(H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad)){
    errorf("Error reading data %s\n",label);
    return(NULL);
  }
  
  /*close dataset*/
  if(H5Dclose(dset)){
    errorf("Error closing data %s\n",label);
    return(NULL);
  }
  if(H5Sclose(space)){
    errorf("Error closing space %s\n",label);
    return(NULL);
  }
  if(H5Tclose(filetype)){
    errorf("Error closing type %s\n",label);
    return(NULL);
  }
  
  TIDY(dims);
  return(jimlad);
}/*read15dCharHDF5*/


/*####################################*/
/*read 1D uint8 array from HDF5*/

uint8_t *read1dUint8HDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  uint8_t *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  //if((filetype!=H5T_NATIVE_USHORT)&&(filetype!=H5T_STD_I16BE)&&(filetype!=H5T_STD_I16LE)){
  //  errorf("Wrong data type\n");
  //  exit(1);
  //}
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return NULL;
  }
  *nBins=dims[0];
  if(!(jimlad=(uint8_t *)calloc(*nBins,sizeof(uint8_t)))){
    errorf("error in float buffer allocation.\n");
    return NULL;
  }
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return NULL;
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  return(jimlad);
}/*read1dUint8HDF5*/


/*####################################*/
/*read 1D uint16 array from HDF5*/

uint16_t *read1dUint16HDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  uint16_t *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  //if((filetype!=H5T_NATIVE_USHORT)&&(filetype!=H5T_STD_I16BE)&&(filetype!=H5T_STD_I16LE)){
  //  errorf("Wrong data type\n");
  //  exit(1);
  //}
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  if(!(jimlad=(uint16_t *)calloc(*nBins,sizeof(uint16_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dUint16HDF5*/


/*####################################*/
/*read 1D uint32 array from HDF5*/

uint32_t *read1dUint32HDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  uint32_t *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  if(!(jimlad=(uint32_t *)calloc(*nBins,sizeof(uint32_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dUint32HDF5*/


/*####################################*/
/*read 1D uint64 array from HDF5*/

uint64_t *read1dUint64HDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  uint64_t *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  if(!(jimlad=(uint64_t *)calloc(*nBins,sizeof(uint64_t)))){
    errorf("error in float buffer allocation.\n");
    return(NULL);
  }
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dUint64HDF5*/


/*#####################################*/
/*read 1D int array from HDF5*/

int *read1dIntHDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  int *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  ASSIGN_CHECKNULL_RETNULL(jimlad,ialloc(dims[0],"",0));
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dIntHDF5*/


/*#####################################*/
/*read 1D float array from HDF5*/

float *read1dFloatHDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  float *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);

  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  ASSIGN_CHECKNULL_RETNULL(jimlad,falloc((uint64_t)dims[0],varName,0));
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }

  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dFloatHDF5*/


/*#####################################*/
/*read 1D double array from HDF5*/

double *read1dDoubleHDF5(hid_t file,char *varName,int *nBins)
{
  int ndims=0;
  hid_t dset,space,filetype;         /* Handles */
  herr_t status;
  hsize_t dims[1];
  double *jimlad=NULL;

  dset=H5Dopen2(file,varName,H5P_DEFAULT);
  filetype=H5Dget_type(dset);
  space=H5Dget_space(dset);
  ndims=H5Sget_simple_extent_dims(space,dims,NULL);
  if(ndims>1){
    errorf("Wrong number of dimensions %d\n",ndims);
    return(NULL);
  }
  *nBins=dims[0];
  ASSIGN_CHECKNULL_RETNULL(jimlad,dalloc(dims[0],"",0));
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    errorf("Data reading error %d\n",status);
    return(NULL);
  }
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if(H5Tclose(filetype)){
    errorf("Error closing type\n");
    return(NULL);
  }
  return(jimlad);
}/*read1dDoubleHDF5*/


/*####################################################*/
/*write a 1D uint8 array*/

int writeComp1dUint8HDF5(hid_t file,char *varName,uint8_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UCHAR);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dUint8HDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

int writeComp1dUint32HDF5(hid_t file,char *varName,uint32_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UINT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dUint32HDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

int writeComp1dInt8HDF5(hid_t file,char *varName,int8_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_CHAR);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dInt8HDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

int writeComp1dInt32HDF5(hid_t file,char *varName,int32_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_INT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dInt32HDF5*/


/*####################################################*/
/*write a 1D uint64 array*/

int writeComp1dUint64HDF5(hid_t file,char *varName,uint64_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_ULONG);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dUint64HDF5*/


/*####################################################*/
/*write a 1D uint16 array*/

int writeComp1dUint16HDF5(hid_t file,char *varName,uint16_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_USHORT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dUint16HDF5*/


/*####################################################*/
/*write a 1D uint32 array*/

int write1dUint32HDF5(hid_t file,char *varName,uint32_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_UINT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);

  return(0);
}/*write1dUint16HDF5*/


/*####################################################*/
/*write a 1D double array*/

int write1dDoubleHDF5(hid_t file,char *varName,double *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_DOUBLE);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write1dDoubleHDF5*/


/*####################################################*/
/*write a 1D char array*/

int write2dCharHDF5(hid_t file,char *varName,char *data,int nWaves,int nBins)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[2];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dims[1]=(hsize_t)nBins;
  dataspace=H5Screate_simple(2,dims,NULL);
  /*datatype=H5Tcopy(H5T_NATIVE_CHAR);*/
  datatype=H5Tcopy(H5T_C_S1);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write2dCharHDF5*/


/*####################################################*/
/*write a uint16 float array*/

int write2dUint16HDF5(hid_t file,char *varName,uint16_t *data,int nWaves,int nBins)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[2];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dims[1]=(hsize_t)nBins;
  dataspace=H5Screate_simple(2,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_USHORT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write2dUint16HDF5*/


/*####################################################*/
/*write a 2D float array*/

int write2dFloatHDF5(hid_t file,char *varName,float *data,int nWaves,int nBins)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[2];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dims[1]=(hsize_t)nBins;
  dataspace=H5Screate_simple(2,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_FLOAT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write2dFloatHDF5*/


/*####################################################*/
/*write a compressed 2D int8 array*/

int writeComp2dInt8HDF5(hid_t file,char *varName,int8_t *data,int nWaves,int nBins)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[2];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[2];


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dims[1]=(hsize_t)nBins;
  dataspace=H5Screate_simple(2,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_CHAR);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  chunk[1]=nBins;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,2,chunk);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp2dInt8HDF5*/


/*####################################################*/
/*write a compressed 2D float array*/

int writeComp2dFloatHDF5(hid_t file,char *varName,float *data,int nWaves,int nBins)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[2];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[2];


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dims[1]=(hsize_t)nBins;
  dataspace=H5Screate_simple(2,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_FLOAT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  chunk[1]=nBins;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,2,chunk);

  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*writeComp2dFloatHDF5*/


/*####################################################*/
/*write a 1D float array*/

int write1dFloatHDF5(hid_t file,char *varName,float *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_FLOAT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write1dFloatHDF5*/


/*####################################################*/
/*write a compressed 1D float array*/

int writeComp1dFloatHDF5(hid_t file,char *varName,float *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_FLOAT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*writeComp1dFloatHDF5*/


/*####################################################*/
/*write a compressed 1D float array*/

int writeComp1dDoubleHDF5(hid_t file,char *varName,double *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/
  hsize_t chunk[1];

  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_DOUBLE);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);

  /*set compression*/
  chunk[0]=nWaves;
  dcpl_id=H5Pcreate (H5P_DATASET_CREATE);
  status=H5Pset_deflate (dcpl_id, 9);
  status=H5Pset_chunk(dcpl_id,1,chunk);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*writeComp1dDoubleHDF5*/


/*####################################################*/
/*write a 1D int array*/

int write1dIntHDF5(hid_t file,char *varName,int *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_INT);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return(-1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return(0);
}/*write1dIntHDF5*/


/*####################################################*/
/*write a 1D int array*/

int write1dInt64HDF5(hid_t file,char *varName,int64_t *data,int nWaves)
{
  hid_t dset;
  herr_t status;
  hsize_t dims[1];
  hid_t datatype,dataspace;  /*data definitions*/
  hid_t lcpl_id,dcpl_id,dapl_id;     /*creation and access properties*/


  /*define dataspace*/
  dims[0]=(hsize_t)nWaves;
  dataspace=H5Screate_simple(1,dims,NULL);
  datatype=H5Tcopy(H5T_NATIVE_LONG);
  /*access and creation properties*/
  lcpl_id=H5Pcopy(H5P_DEFAULT);
  dcpl_id=H5Pcopy(H5P_DEFAULT);
  dapl_id=H5Pcopy(H5P_DEFAULT);


  /*create new dataset*/
  dset=H5Dcreate2(file,varName,datatype,dataspace,lcpl_id,dcpl_id,dapl_id);
  if(dset<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    errorf("Error writing %s\n",varName);
    return -1;
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return 0;
}/*write1dInt64HDF5*/


/*the end*/
/*####################################################*/

