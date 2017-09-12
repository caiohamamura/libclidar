#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "hdf5.h"
#include "stdint.h"
#include "tools.h"
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

lvisLGWdata *readLVISlgw(char *namen,int *nWaves)
{
  uint64_t i=0,len=0;
  uint64_t offset=0;
  lvisLGWdata *data=NULL;
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
  *nWaves=(int)(len/sizeof(lvisLGWdata));
  if(!(data=(lvisLGWdata *)calloc(*nWaves,sizeof(lvisLGWdata)))){
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
}/*readLVISlgw*/


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
  void checkNumber(int,int,char *);
  float *read1dFloatHDF5(hid_t,char *,int *);
  double *read1dDoubleHDF5(hid_t,char *,int *);
  uint32_t *read1dUint32HDF5(hid_t,char *,int *);
  uint16_t **readDatasetHDF5(hid_t,char *,int *,int *);
  hid_t file;         /* Handles */


  /*allocate structure*/
  if(!(lvis=(lvisHDF *)calloc(1,sizeof(lvisHDF)))){
    fprintf(stderr,"error in LVIS structure allocation.\n");
    exit(1);
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
  fprintf(stdout,"Reading %s\n",inNamen);
  file=H5Fopen(inNamen,H5F_ACC_RDONLY,H5P_DEFAULT);

  /*read 1D double arrays*/
  lvis->lon0=read1dDoubleHDF5(file,"LON0",&nWaves);
  lvis->nWaves=nWaves;
  lvis->lat0=read1dDoubleHDF5(file,"LAT0",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"LAT0");
  lvis->lon1023=read1dDoubleHDF5(file,"LON1023",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"LON1023");
  lvis->lat1023=read1dDoubleHDF5(file,"LAT1023",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"LAT1023");
  lvis->time=read1dDoubleHDF5(file,"TIME",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"TIME");

  /*read 1D float arrays*/
  lvis->zen=read1dFloatHDF5(file,"INCIDENTANGLE",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"INCIDENTANGLE");
  lvis->z0=read1dFloatHDF5(file,"Z0",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"Z0");
  lvis->z1023=read1dFloatHDF5(file,"Z1023",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"Z1023");
  lvis->sigmean=read1dFloatHDF5(file,"SIGMEAN",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"SIGMEAN");

  /*read 1D uint32 arrays*/
  lvis->lfid=read1dUint32HDF5(file,"LFID",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"LFID");
  lvis->shotN=read1dUint32HDF5(file,"SHOTNUMBER",&nWaves);
  checkNumber(nWaves,lvis->nWaves,"SHOTNUMBER");

  /*read 2d unit16 arrays*/
  lvis->wave=readDatasetHDF5(file,"RXWAVE",&lvis->nBins,&nWaves);
  checkNumber(nWaves,lvis->nWaves,"RXWAVE");
  lvis->pulse=readDatasetHDF5(file,"TXWAVE",&lvis->pBins,&nWaves);
  checkNumber(nWaves,lvis->nWaves,"TXWAVE");

  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }

  return(lvis);
}/*readLVIShdf*/


/*#####################################*/
/*check integers match*/

void checkNumber(int newNumb,int oldNumb,char *label)
{
  if(newNumb!=oldNumb){
    fprintf(stderr,"Number mismatch %d %d for %s\n",newNumb,oldNumb,label);
    exit(1);
  }
  return;
}/*checkNumber*/


/*#####################################*/
/*read a HDF5 dataset*/

uint16_t **readDatasetHDF5(hid_t file,char *label,int *nBins,int *nWaves)
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
    fprintf(stderr,"error in float buffer allocation.\n");
    exit(1);
  }

  if(H5Sget_simple_extent_dims(space,dims,NULL)!=ndims){
    fprintf(stderr,"Error\n");
    exit(1);
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
    fprintf(stderr,"Error reading data %s\n",label);
    exit(1);
  }

  /*close dataset*/
  if(H5Dclose(dset)){
    fprintf(stderr,"Error closing data %s\n",label);
    exit(1);
  }
  if(H5Sclose(space)){
    fprintf(stderr,"Error closing space %s\n",label);
    exit(1);
  }

  TIDY(dims);
  return(jimlad);
}/*readDatasetHDF5*/


/*#####################################*/
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
    fprintf(stderr,"Wrong number of dimensions %d\n",ndims);
    exit(1);
  }
  *nBins=dims[0];
  if(!(jimlad=(uint32_t *)calloc(*nBins,sizeof(uint32_t)))){
    fprintf(stderr,"error in float buffer allocation.\n");
    exit(1);
  }
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    fprintf(stderr,"Data reading error %d\n",status);
    exit(1);
  }
  status=H5Dclose(dset);

  return(jimlad);
}/*read1dUint32HDF5*/


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
    fprintf(stderr,"Wrong number of dimensions %d\n",ndims);
    exit(1);
  }
  *nBins=dims[0];
  jimlad=falloc(dims[0],"",0);
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    fprintf(stderr,"Data reading error %d\n",status);
    exit(1);
  }
  status=H5Dclose(dset);

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
    fprintf(stderr,"Wrong number of dimensions %d\n",ndims);
    exit(1);
  }
  *nBins=dims[0];
  jimlad=dalloc(dims[0],"",0);
  status=H5Dread(dset,filetype,H5S_ALL,H5S_ALL,H5P_DEFAULT,jimlad);
  if(status){
    fprintf(stderr,"Data reading error %d\n",status);
    exit(1);
  }
  status=H5Dclose(dset);

  return(jimlad);
}/*read1dDoubleHDF5*/


/*####################################################*/
/*write a 1D float array*/

void write1dDoubleHDF5(hid_t file,char *varName,double *data,int nWaves)
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
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  return;
}/*write1dDoubleHDF5*/


/*####################################################*/
/*write a 1D char array*/

void write2dCharHDF5(hid_t file,char *varName,char *data,int nWaves,int nBins)
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
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  return;
}/*write2dCharHDF5*/


/*####################################################*/
/*write a 2D float array*/

void write2dFloatHDF5(hid_t file,char *varName,float *data,int nWaves,int nBins)
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
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  status=H5Sclose(dataspace);
  return;
}/*write2dFloatHDF5*/


/*####################################################*/
/*write a 1D float array*/

void write1dFloatHDF5(hid_t file,char *varName,float *data,int nWaves)
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
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  return;
}/*write1dFloatHDF5*/


/*####################################################*/
/*write a 1D float array*/

void write1dIntHDF5(hid_t file,char *varName,int *data,int nWaves)
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
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*write data*/
  status=H5Dwrite(dset,datatype,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void *)data);
  if(status<0){
    fprintf(stderr,"Error writing %s\n",varName);
    exit(1);
  }

  /*close data*/
  status=H5Dclose(dset);
  return;
}/*write1dIntHDF5*/


/*####################################################*/
/*write data to HDF5*/

void writeGEDIhdf(gediHDF *hdfData,char *namen)
{ 
  hid_t file;         /* Handles */
  
  /*open new file*/
  file=H5Fcreate(namen,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  
  /*write header*/
  write1dIntHDF5(file,"NWAVES",&hdfData->nWaves,1);
  write1dIntHDF5(file,"NBINS",&hdfData->nBins,1);
  write1dFloatHDF5(file,"PSIGMA",&hdfData->pSigma,1);
  write1dFloatHDF5(file,"FSIGMA",&hdfData->fSigma,1);
  /*write datasets*/
  write1dDoubleHDF5(file,"LON0",hdfData->lon,hdfData->nWaves);
  write1dDoubleHDF5(file,"LAT0",hdfData->lat,hdfData->nWaves);
  write1dFloatHDF5(file,"SLOPE",hdfData->slope,hdfData->nWaves);
  write1dFloatHDF5(file,"COVER",hdfData->cov,hdfData->nWaves);
  write1dFloatHDF5(file,"ZG",hdfData->gElev,hdfData->nWaves);
  write1dFloatHDF5(file,"ZGdem",hdfData->demElev,hdfData->nWaves);
  write1dFloatHDF5(file,"BEAMDENSE",hdfData->beamDense,hdfData->nWaves);
  write1dFloatHDF5(file,"POINTDENSE",hdfData->pointDense,hdfData->nWaves);
  write1dFloatHDF5(file,"INCIDENTANGLE",hdfData->zen,hdfData->nWaves);
  write2dFloatHDF5(file,"RXWAVE",hdfData->wave,hdfData->nWaves,hdfData->nBins);
  if(hdfData->ground)write2dFloatHDF5(file,"GRWAVE",hdfData->ground,hdfData->nWaves,hdfData->nBins);
  write1dFloatHDF5(file,"Z0",hdfData->z0,hdfData->nWaves);
  write1dFloatHDF5(file,"ZN",hdfData->zN,hdfData->nWaves);
  write2dCharHDF5(file,"WAVEID",hdfData->waveID,hdfData->nWaves,hdfData->idLength);
  
  /*close file*/
  if(H5Fclose(file)){
    fprintf(stderr,"Issue closing file\n");
    exit(1);
  }
  fprintf(stdout,"Waveforms written to %s\n",namen);
  return;
}/*writeGEDIhdf*/


/*####################################################*/
/*tidy GEDI HDF data structire*/

gediHDF *tidyGediHDF(gediHDF *hdfData)
{

  if(hdfData){
    TIDY(hdfData->wave);
    TIDY(hdfData->ground);
    TIDY(hdfData->waveID);
    TIDY(hdfData->z0);       /*wave top elevations*/
    TIDY(hdfData->zN);       /*wave bottom elevations*/
    TIDY(hdfData->lon);     /*longitudes*/
    TIDY(hdfData->lat);     /*latitudes*/
    TIDY(hdfData->slope);    /*ground slope*/
    TIDY(hdfData->cov);      /*canopy cover*/
    TIDY(hdfData->gElev);    /*ground elevation, CofG*/
    TIDY(hdfData->demElev);  /*ground elevation, DEM*/
    TIDY(hdfData->beamDense);/*beam density*/
    TIDY(hdfData->pointDense);/*point density*/
    TIDY(hdfData->zen);      /*scan angles, or mean angles*/
    TIDY(hdfData);
  }

  return(hdfData);
}/*tidyHDFdata*/


/*the end*/
/*####################################################*/

