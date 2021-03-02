#include "Temperature/FindQTemp.h"

////////////////////////////////////////////////////////////////////////////////
void findTempPress(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0};
   int   ey[6]={   0,  0,  1, -1,  0,  0};
   int   ez[6]={   0,  0,  0,  0,  1, -1};
   real T, test;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   int nx                                = para->getParH(para->getCoarse())->nx;
   int ny                                = para->getParH(para->getCoarse())->ny;
   unsigned int nnx                      = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny                      = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz                      = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                          = para->getParH(para->getCoarse())->geo;
   unsigned int* kk                      = para->getParH(para->getCoarse())->k;
   // real TempBC                        = para->getTemperatureBC();
   // real VelBC                         = para->getVelocity();
   TempPressforBoundaryConditions  Temp  = para->getParH(para->getCoarse())->TempPress;
   Temp.kTemp = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   k = nnz + STARTOFFZ - 2;
   {
   //for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){          
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test = (real)-1.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     if (ez[l]==1)
                     {
                        T = (real)1.f;
                     } 
                     //else if (ez[l]==-1)
                     //{
                     //   T = (real)1.f;//2.f;
                     //}
                     test = (real)1.f;
                  }
               }
               if (test == (real)1.f)
               {
                  if (T == (real)1.f)
                  {
                     Temp.k[Temp.kTemp]          = kk[m];
                     Temp.temp[Temp.kTemp]       = (real)0.f;//TempBC;                  
                     Temp.velo[Temp.kTemp]       = (real)0.f;//VelBC;                  
                     Temp.kTemp++;              
                  }                         
                  //else if (T == (real)2.f)
                  //{
                  //   Temp.k[Temp.kTemp]          = kk[m];
                  //   Temp.temp[Temp.kTemp]       = -TempBC;                  
                  //   Temp.velo[Temp.kTemp]       = -VelBC;                  
                  //   Temp.kTemp++;              
                  //}
                  //else
                  //{
                  //   Temp.k[Temp.kTemp]          = kk[m];
                  //   Temp.temp[Temp.kTemp]       = 0.0;                  
                  //   Temp.kTemp++;              
                  //}
               }
            }     
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void findKforTempPress(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B   ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0,};
   int   ey[6]={   0,  0,  1, -1,  0,  0,};
   int   ez[6]={   0,  0,  0,  0,  1, -1,};
   real ON;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   real test = (real)0.f;
   int nx                                = para->getParH(para->getCoarse())->nx;
   int ny                                = para->getParH(para->getCoarse())->ny;
   unsigned int nnx                      = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny                      = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz                      = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                          = para->getParH(para->getCoarse())->geo;
   TempPressforBoundaryConditions  Temp  = para->getParH(para->getCoarse())->TempPress;
   Temp.kTemp = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   k = nnz + STARTOFFZ - 2;
   {
   //for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){       
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test =(real)0.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     ON =(real) 1.f; 
                  }
                  else{
                     ON = (real)0.f;
                  }
                  test += ON;
               }
               if (test>0)
               {
                  Temp.kTemp++;                         
               }                         
            }     
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
void findTempVel(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0};
   int   ey[6]={   0,  0,  1, -1,  0,  0};
   int   ez[6]={   0,  0,  0,  0,  1, -1};
   real T, test;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   int nx                             = para->getParH(para->getCoarse())->nx;
   int ny                             = para->getParH(para->getCoarse())->ny;
   unsigned int nnx                   = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny                   = para->getParH(para->getCoarse())->gridNY;
   // unsigned int nnz                   = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                       = para->getParH(para->getCoarse())->geo;
   unsigned int* kk                   = para->getParH(para->getCoarse())->k;
   real TempBC                     = para->getTemperatureBC();
   real VelBC                      = para->getVelocity();
   TempVelforBoundaryConditions  Temp = para->getParH(para->getCoarse())->TempVel;
   Temp.kTemp = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   k = STARTOFFZ + 1;
   {
      //for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){          
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test = (real)-1.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     /*if (ez[l]==1)
                     {
                        T = (real)1.f;
                     } 
                     else*/ if (ez[l]==-1)
                     {
                        T = (real)1.f;//2.f;
                     }
                     test = (real)1.f;
                  }
               }
               if (test == (real)1.f)
               {
                  if (T == (real)1.f)
                  {
                     Temp.k[Temp.kTemp]          = kk[m];
                     Temp.temp[Temp.kTemp]       = TempBC;                  
                     Temp.velo[Temp.kTemp]       = VelBC;                  
                     Temp.kTemp++;              
                  }                         
                  //else if (T == (real)2.f)
                  //{
                  //   Temp.k[Temp.kTemp]          = kk[m];
                  //   Temp.temp[Temp.kTemp]       = -TempBC;                  
                  //   Temp.velo[Temp.kTemp]       = -VelBC;                  
                  //   Temp.kTemp++;              
                  //}
                  //else
                  //{
                  //   Temp.k[Temp.kTemp]          = kk[m];
                  //   Temp.temp[Temp.kTemp]       = 0.0;                  
                  //   Temp.kTemp++;              
                  //}
               }
            }     
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void findKforTempVel(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B   ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0,};
   int   ey[6]={   0,  0,  1, -1,  0,  0,};
   int   ez[6]={   0,  0,  0,  0,  1, -1,};
   real ON;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   real test = (real)0.f;
   int nx                             = para->getParH(para->getCoarse())->nx;
   int ny                             = para->getParH(para->getCoarse())->ny;
   unsigned int nnx                   = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny                   = para->getParH(para->getCoarse())->gridNY;
   // unsigned int nnz                   = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                       = para->getParH(para->getCoarse())->geo;
   TempVelforBoundaryConditions Temp  = para->getParH(para->getCoarse())->TempVel;
   Temp.kTemp = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   k = STARTOFFZ + 1;
   {
      //for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){       
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test =(real)0.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     ON =(real) 1.f; 
                  }
                  else{
                     ON = (real)0.f;
                  }
                  test += ON;
               }
               if (test>0)
               {
                  Temp.kTemp++;                         
               }                         
            }     
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void findTemp(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0};
   int   ey[6]={   0,  0,  1, -1,  0,  0};
   int   ez[6]={   0,  0,  0,  0,  1, -1};
   real ON[7];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   int nx                          = para->getParH(para->getCoarse())->nx;
   int ny                          = para->getParH(para->getCoarse())->ny;
   int *geo_mat                    = para->getParH(para->getCoarse())->geo;
   unsigned int nnx                = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny                = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz                = para->getParH(para->getCoarse())->gridNZ;
   unsigned int* kk                = para->getParH(para->getCoarse())->k;
   real TempBC                  = para->getTemperatureBC();
   TempforBoundaryConditions Temp = para->getParH(para->getCoarse())->Temp;
   Temp.kTemp = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){          
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               ON[6] = (real)-1.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     if (ey[l]==1)
                     {
                        ON[2] = (real)1.f;
                     } 
                     else if (ey[l]==-1)
                     {
                        ON[2] = (real)2.f;
                     }
                     ON[6] = (real)1.f;
                  }
               }
               if (ON[6] == (real)1.f)
               {
                  if (ON[2] == (real)1.f)
                  {
                     Temp.k[Temp.kTemp]          = kk[m];
                     Temp.temp[Temp.kTemp]       = TempBC;                  
                     Temp.kTemp++;              
                  }                         
                  else if (ON[2] == (real)2.f)
                  {
                     Temp.k[Temp.kTemp]          = kk[m];
                     Temp.temp[Temp.kTemp]       = -TempBC;                  
                     Temp.kTemp++;              
                  }
                  else
                  {
                     Temp.k[Temp.kTemp]          = kk[m];
                     Temp.temp[Temp.kTemp]       = 0.0;                  
                     Temp.kTemp++;              
                  }
               }
            }     
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void findKforTemp(Parameter* para)
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B   ////////////////////////
   int   ex[6]={   1, -1,  0,  0,  0,  0,};
   int   ey[6]={   0,  0,  1, -1,  0,  0,};
   int   ez[6]={   0,  0,  0,  0,  1, -1,};
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   real test = (real)0.f;
   real ON;
   para->getTempH()->kTemp = 0;
   int nx           = para->getParH(para->getCoarse())->nx;
   int ny           = para->getParH(para->getCoarse())->ny;
   unsigned int nnx = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz = para->getParH(para->getCoarse())->gridNZ;
   int *geo_mat     = para->getParH(para->getCoarse())->geo;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(k=STARTOFFZ+1; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY+1; j<=nny+STARTOFFY-2 ; j++){       
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test =(real)0.f;
               for(l=0;l<=5;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     ON =(real) 1.f; 
                  }
                  else{
                     ON = (real)0.f;
                  }
                  test += ON;
               }
               if (test>0)
               {
                  para->getParH(para->getCoarse())->Temp.kTemp++;                         
               }                         
            }     
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

