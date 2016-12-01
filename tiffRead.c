#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "geotiffio.h"
#include "xtiffio.h"
#include "tools.h"
#include "tiffRead.h"

/*############################################*/
/*read a geotiff*/

void readGeotiff(geot *geotiff,char *namen,char readData)
{
  int i=0;
  unsigned int type=0;
  short tiepointsize=0,scalesize=0;
  double *tiepoints=NULL,*scale=NULL;
  //GTIF *tiffStruct=(GTIF*)0;
  TIFF *tiffIn=(TIFF*)0;

  if((tiffIn=XTIFFOpen(namen,"r"))==NULL){
    fprintf(stderr,"Damn, no %s\n",namen);
    exit(1);
  }
  geotiff->tiepoints=dalloc(6,"",0);
  geotiff->scale=dalloc(3,"",0);

  /*for now we have assumed OSNG, 27700, projection*/
  //tiffStruct=GTIFNew(tiffIn);
  TIFFGetField(tiffIn,TIFFTAG_IMAGEWIDTH,&geotiff->nX);
  TIFFGetField(tiffIn,TIFFTAG_IMAGELENGTH,&geotiff->nY);
  TIFFGetField(tiffIn,TIFFTAG_DATATYPE,&type);

  //GTIFPrint(tiffStruct,0,0);

  TIFFGetField(tiffIn,TIFFTAG_GEOPIXELSCALE,&scalesize,&scale);
  TIFFGetField(tiffIn,TIFFTAG_GEOTIEPOINTS,&tiepointsize,&tiepoints);
  for(i=0;i<3;i++) geotiff->scale[i]=scale[i];
  for(i=0;i<6;i++) geotiff->tiepoints[i]=tiepoints[i];

  if(readData){
    if(type==0){  /*unsigned char*/
      geotiff->image=uchalloc(geotiff->nX*geotiff->nY,namen,0);
      for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
        if(TIFFReadScanline(tiffIn,&(geotiff->image[i*geotiff->nX]),i,1)!=1){
          fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
          exit(1);
        }
      }
    }else if(type==3){ /*float*/
      if(((int)TIFFScanlineSize(tiffIn)/geotiff->nX)==4){
        geotiff->fImage=falloc(geotiff->nX*geotiff->nY,namen,0);
        for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
          if(TIFFReadScanline(tiffIn,&(geotiff->fImage[i*geotiff->nX]),i,1)!=1){
            fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
            exit(1);
          }
        }
      }else if(((int)TIFFScanlineSize(tiffIn)/geotiff->nX)==8){
        geotiff->dImage=dalloc(geotiff->nX*geotiff->nY,namen,0);
        for(i=0;i<geotiff->nY;i++){                  /*looping along the lattitude*/
          if(TIFFReadScanline(tiffIn,&(geotiff->dImage[i*geotiff->nX]),i,1)!=1){
            fprintf(stderr,"Error reading scan line %d from tiff image\n",i);
            exit(1);
          }
        }
      }else{
        fprintf(stderr,"What do you think you're doing!?!\n");
        exit(1);
      }
    }else{
      fprintf(stderr,"Cannot handle type %d\n",type);
      exit(1);
    }
  }/*read data question*/

  XTIFFClose(tiffIn);
  return;
}/*readGeotiff*/


/*############################################*/
/*tidy tiff file structure*/

geot *tidyTiff(geot *tiff)
{
  if(tiff){
    TIDY(tiff->image);
    TIDY(tiff->fImage);
    TIDY(tiff->dImage);
    TIDY(tiff->tiepoints);
    TIDY(tiff->scale);
    TIDY(tiff);
  }

  return(tiff);
}/*tidyTiff*/
/*the end*/
/*############################################*/

