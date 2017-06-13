

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
/*LVIS v1.03 structure*/

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
}lvisData103;


/*#######################################*/
/*functions*/

lvisData103 *readData(char *,int *);
void checkLVISsizes();

/*the end*/
/*####################################################*/

