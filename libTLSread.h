
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


/*##########################################*/
/*TLS data*/

typedef struct{
  int bin;       /*bin number*/
  float x;       /*coordinate*/
  float y;       /*coordinate*/
  float z;       /*coordinate*/
  float gap;     /*voxel gap fraction*/
  float r;       /*range*/
  uint16_t refl; /*reflectance*/
}tlsPoint;

/*the end*/
/*##########################################*/

