/*########################*/
/*# structures for voxel #*/
/*# lidar programs       #*/
/*########################*/


/*###################################################*/
/*voxels*/

typedef struct{
  int nScans;        /*number of contributing files*/
  float **hits;      /*hits per scan location*/
  float **miss;      /*misses per scan location*/
  float *rmse;       /*rmse of signal going in*/
  int *contN;        /*for normalising ALS*/
  int nVox;          /*total number of voxels*/
  int nX;            /*number of voxels*/
  int nY;            /*number of voxels*/
  int nZ;            /*number of voxels*/
  double res[3];     /*voxel resolution*/
  double bounds[6];  /*voxel bounds minX minY minZ maxX maxY maxZ*/
  char useRMSE;      /*switch to save RAM in voxelate*/
}voxStruct;


/*###################################################*/
/*lidar radiometric parameters for voxels*/

typedef struct{
  float minRefl;   /*minimum refletance value to scale between 0 and 1*/
  float maxRefl;   /*maximum refletance value to scale between 0 and 1*/
  float appRefl;   /*scale between TLS reflectance and size*/
  float beamTanDiv; /*tan of tls beam divergence*/
  float beamRad;    /*TLS start radius*/
  float minGap;    /*minimum gap fraction correction to apply*/
}lidVoxPar;


/*#####################################*/
/*structure to hold a range image*/

typedef struct{
  int nBins;        /*number of range bins*/
  int nX;           /*number of x bins*/
  int nY;           /*number of y bins*/
  float rRes;       /*range resolution*/
  float iRes;       /*image resolution*/
  double bounds[6]; /*minX minY minZ, maxX maxY maxZ*/
  double x0;        /*central coordinate*/
  double y0;        /*central coordinate*/
  double z0;        /*central coordinate*/
  char **image;     /*gap image at each range*/
  float grad[3];    /*vector along rage image*/
}rImageStruct;

/*the end*/
/*###################################################*/

