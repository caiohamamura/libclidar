

/*#######################*/
/*# A library for       #*/
/*# handling LVIS files #*/
/*# S Hancock, 2017     #*/
/*#######################*/

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




/*#######################################*/
/*LVIS LGW v1.03 and below structure*/

typedef struct{
   uint32_t lfid;       /* LVIS file identifier*/
   uint32_t shotN; /* LVIS shotnumber*/
   float az;            /* true heading from the aircraft to the ground (degrees)*/
   float zen;      /*zenith angle (deghrees)*/
   float range;      /* range from the aircraft to the ground (meters)*/
   double lvistime;  /* LVIS recorded UTC time (seconds of the day) when the shot was acquired*/
   double lon0;      /* longitude of the highest sample of the waveform (degrees east)*/
   double lat0;      /* latitude of the highest sample of the waveform (degrees north)*/
   float  z0;        /* elevation of the highest sample of the waveform (m)*/
   double lon431;    /* longitude of the lowest sample of the waveform (degrees east)*/
   double lat431;    /* latitude of the lowest sample of the waveform (degrees north)*/
   float  z431;      /* elevation of the lowest sample of the waveform (m)*/
   float  sigmean;   /* signal mean noise level, calculated in-flight (counts)*/
   unsigned char txwave[80];  /* transmit waveform, recorded in-flight (counts)*/
   unsigned char rxwave[432]; /* return   waveform, recorded in-flight (counts)*/
}lvisLGWdata;


/*#######################################*/
/*LVIS overall structure*/

typedef struct{
  int verMaj;      /*major version*/
  int verMin;      /*minor version*/
  int nWaves;      /*number of waveforms*/
  int nBins;       /*number of waveform bins*/
  FILE *ipoo;      /*input file*/
  lvisLGWdata *data;  /*data pointer*/
  char byteord;    /*byte order of this computer*/
}lvisLGWstruct;


/*#####################################*/
/*LVIS HDF5 structure*/

typedef struct{
  int nWaves;   /*number of waveforms*/
  int nBins;    /*number of waveform bins*/
  int pBins;     /*number of pulse bins*/
  /*data per wave*/
  double *lon0;       /*LON0*/
  double *lat0;       /*LAT0*/
  double *lon1023;    /*LON1023*/
  double *lat1023;    /*LAT1023*/
  uint32_t *lfid;     /*LFID*/
  uint32_t *shotN;    /*SHOTNUMBER*/
  uint16_t **wave;    /*RXWAVE*/
  uint16_t **pulse;   /*TXWAVE*/
  float *zen;         /*INCIDENTANGLE*/
  float *z0;         /*Z0*/
  float *z1023;      /*Z1023*/
  float *sigmean;     /*SIGMEAN*/
  double *time;       /*TIME*/
}lvisHDF;


/*#######################################*/
/*functions*/

lvisLGWdata *readLVISlgw(char *,int *);
void checkLVISsizes();
lvisHDF *tidyLVISstruct(lvisHDF *);
lvisHDF *readLVIShdf(char *);

/*the end*/
/*####################################################*/

