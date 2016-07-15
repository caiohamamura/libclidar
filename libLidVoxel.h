/*########################*/
/*# structures for voxel #*/
/*# lidar programs       #*/
/*########################*/


/*#####################################*/
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
  double bounds[6];  /*voxel bounds minX maxX minY maxY minZ maxZ*/
  char useRMSE;      /*switch to save RAM in voxelate*/
}voxStruct;

/*the end*/
/*#####################################*/

