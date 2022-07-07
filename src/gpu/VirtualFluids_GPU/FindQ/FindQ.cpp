#include "FindQ/FindQ.h"
#include "logger/Logger.h"
#include "lbm/constants/D3Q27.h"

using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
void findQ(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

   VF_LOG_CRITICAL("findQ() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   int nx                       = para->getParH(lev)->nx;
   int ny                       = para->getParH(lev)->ny;
   unsigned int nnx             = para->getParH(lev)->gridNX;
   unsigned int nny             = para->getParH(lev)->gridNY;
   unsigned int nnz             = para->getParH(lev)->gridNZ;
   int* geo_mat                 = para->getParH(lev)->geo;
   unsigned int* kk             = para->getParH(para->getCoarse())->k;
   unsigned int sizeQ           = para->getParH(lev)->noSlipBC.numberOfBCnodes;
   real* QQ                     = para->getParH(lev)->noSlipBC.q27[0];
   QforBoundaryConditions &QIN  = para->getParH(lev)->noSlipBC;
   QIN.numberOfBCnodes = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //int relx, rely, relz;

   //unsigned int centerX = nnx / 2;
   //unsigned int centerY = nny / 2;
   //unsigned int centerZ = nnz / 2;
   //real        radius  = nny / 5.f;//2.56f;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   QforBoundaryConditions Q;
   Q.q27[E   ] = &QQ[E   *sizeQ];
   Q.q27[W   ] = &QQ[W   *sizeQ];
   Q.q27[N   ] = &QQ[N   *sizeQ];
   Q.q27[S   ] = &QQ[S   *sizeQ];
   Q.q27[T   ] = &QQ[T   *sizeQ];
   Q.q27[B   ] = &QQ[B   *sizeQ];
   Q.q27[NE  ] = &QQ[NE  *sizeQ];
   Q.q27[SW  ] = &QQ[SW  *sizeQ];
   Q.q27[SE  ] = &QQ[SE  *sizeQ];
   Q.q27[NW  ] = &QQ[NW  *sizeQ];
   Q.q27[TE  ] = &QQ[TE  *sizeQ];
   Q.q27[BW  ] = &QQ[BW  *sizeQ];
   Q.q27[BE  ] = &QQ[BE  *sizeQ];
   Q.q27[TW  ] = &QQ[TW  *sizeQ];
   Q.q27[TN  ] = &QQ[TN  *sizeQ];
   Q.q27[BS  ] = &QQ[BS  *sizeQ];
   Q.q27[BN  ] = &QQ[BN  *sizeQ];
   Q.q27[TS  ] = &QQ[TS  *sizeQ];
   Q.q27[REST] = &QQ[REST*sizeQ];
   Q.q27[TNE ] = &QQ[TNE *sizeQ];
   Q.q27[TSW ] = &QQ[TSW *sizeQ];
   Q.q27[TSE ] = &QQ[TSE *sizeQ];
   Q.q27[TNW ] = &QQ[TNW *sizeQ];
   Q.q27[BNE ] = &QQ[BNE *sizeQ];
   Q.q27[BSW ] = &QQ[BSW *sizeQ];
   Q.q27[BSE ] = &QQ[BSE *sizeQ];
   Q.q27[BNW ] = &QQ[BNW *sizeQ];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(k=STARTOFFZ + 1 ; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY + 1 ; j<=nny+STARTOFFY-2 ; j++){          //j<=nny/2+STARTOFFY     //j<=STARTOFFY+1
         for(i=STARTOFFX + 1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               //relx = i - STARTOFFX - centerX;
               //rely = j - STARTOFFY - centerY;
               //relz = k - STARTOFFZ - centerZ;
               ON[18] = (real)-1.f;
               for(l=0;l<=26;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     //ON[l] = -(((real)ex[l]*(real)relx + (real)ey[l]*(real)rely + (real)ez[l]*(real)relz +/*+/- Achtung, innen und aussen nicht verwechseln!!*/
                     //           sqrt(pow((real)ex[l]*(real)relx + (real)ey[l]*(real)rely + (real)ez[l]*(real)relz,2) +
                     //               (pow((real)ex[l],2) + pow((real)ey[l],2) + pow((real)ez[l],2))* (pow(radius,2) -
                     //                pow((real)relx,2) - pow((real)rely,2) - pow((real)relz,2))))
                     //           /(pow((real)ex[l],2) + pow((real)ey[l],2) + pow((real)ez[l],2)));
                     ON[l] = (real)0.5f;//1.0f;
                     ON[18] = (real)1.f; //ZERO
                  }
                  else{
                     ON[l] = (real)-1.f;
                  }
               }
               if (ON[18]==(real)1.f)
               {
                  QIN.k[QIN.numberOfBCnodes]          = kk[m];

                  //Q.q27[E   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[W   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[N   ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[S   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[T   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[B   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[NE  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[SW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[SE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[NW  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[TE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[BW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[BE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[TW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[TN  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[BS  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[BN  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[TS  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[REST][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[TNE ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[TSW ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[TSE ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[TNW ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[BNE ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[BSW ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[BSE ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[BNW ][QIN.numberOfBCnodes] = 0.f;

                  //Q.q27[E   ][QIN.numberOfBCnodes] = ON[W   ];
                  //Q.q27[W   ][QIN.numberOfBCnodes] = ON[E   ];
                  //Q.q27[N   ][QIN.numberOfBCnodes] = ON[S   ];
                  //Q.q27[S   ][QIN.numberOfBCnodes] = ON[N   ];
                  //Q.q27[T   ][QIN.numberOfBCnodes] = ON[B   ];
                  //Q.q27[B   ][QIN.numberOfBCnodes] = ON[T   ];
                  //Q.q27[NE  ][QIN.numberOfBCnodes] = ON[SW  ];
                  //Q.q27[SW  ][QIN.numberOfBCnodes] = ON[NE  ];
                  //Q.q27[SE  ][QIN.numberOfBCnodes] = ON[NW  ];
                  //Q.q27[NW  ][QIN.numberOfBCnodes] = ON[SE  ];
                  //Q.q27[TE  ][QIN.numberOfBCnodes] = ON[BW  ];
                  //Q.q27[BW  ][QIN.numberOfBCnodes] = ON[TE  ];
                  //Q.q27[BE  ][QIN.numberOfBCnodes] = ON[TW  ];
                  //Q.q27[TW  ][QIN.numberOfBCnodes] = ON[BE  ];
                  //Q.q27[TN  ][QIN.numberOfBCnodes] = ON[BS  ];
                  //Q.q27[BS  ][QIN.numberOfBCnodes] = ON[TN  ];
                  //Q.q27[BN  ][QIN.numberOfBCnodes] = ON[TS  ];
                  //Q.q27[TS  ][QIN.numberOfBCnodes] = ON[BN  ];
                  //Q.q27[REST][QIN.numberOfBCnodes] = ON[REST];
                  //Q.q27[TNE ][QIN.numberOfBCnodes] = ON[BSW ];
                  //Q.q27[TSW ][QIN.numberOfBCnodes] = ON[BNE ];
                  //Q.q27[TSE ][QIN.numberOfBCnodes] = ON[BNW ];
                  //Q.q27[TNW ][QIN.numberOfBCnodes] = ON[BSE ];
                  //Q.q27[BNE ][QIN.numberOfBCnodes] = ON[TSW ];
                  //Q.q27[BSW ][QIN.numberOfBCnodes] = ON[TNE ];
                  //Q.q27[BSE ][QIN.numberOfBCnodes] = ON[TNW ];
                  //Q.q27[BNW ][QIN.numberOfBCnodes] = ON[TSE ];

                  Q.q27[E   ][QIN.numberOfBCnodes] = ON[E   ];
                  Q.q27[W   ][QIN.numberOfBCnodes] = ON[W   ];
                  Q.q27[N   ][QIN.numberOfBCnodes] = ON[N   ];
                  Q.q27[S   ][QIN.numberOfBCnodes] = ON[S   ];
                  Q.q27[T   ][QIN.numberOfBCnodes] = ON[T   ];
                  Q.q27[B   ][QIN.numberOfBCnodes] = ON[B   ];
                  Q.q27[NE  ][QIN.numberOfBCnodes] = ON[NE  ];
                  Q.q27[SW  ][QIN.numberOfBCnodes] = ON[SW  ];
                  Q.q27[SE  ][QIN.numberOfBCnodes] = ON[SE  ];
                  Q.q27[NW  ][QIN.numberOfBCnodes] = ON[NW  ];
                  Q.q27[TE  ][QIN.numberOfBCnodes] = ON[TE  ];
                  Q.q27[BW  ][QIN.numberOfBCnodes] = ON[BW  ];
                  Q.q27[BE  ][QIN.numberOfBCnodes] = ON[BE  ];
                  Q.q27[TW  ][QIN.numberOfBCnodes] = ON[TW  ];
                  Q.q27[TN  ][QIN.numberOfBCnodes] = ON[TN  ];
                  Q.q27[BS  ][QIN.numberOfBCnodes] = ON[BS  ];
                  Q.q27[BN  ][QIN.numberOfBCnodes] = ON[BN  ];
                  Q.q27[TS  ][QIN.numberOfBCnodes] = ON[TS  ];
                  Q.q27[REST][QIN.numberOfBCnodes] = ON[REST];
                  Q.q27[TNE ][QIN.numberOfBCnodes] = ON[TNE ];
                  Q.q27[TSW ][QIN.numberOfBCnodes] = ON[TSW ];
                  Q.q27[TSE ][QIN.numberOfBCnodes] = ON[TSE ];
                  Q.q27[TNW ][QIN.numberOfBCnodes] = ON[TNW ];
                  Q.q27[BNE ][QIN.numberOfBCnodes] = ON[BNE ];
                  Q.q27[BSW ][QIN.numberOfBCnodes] = ON[BSW ];
                  Q.q27[BSE ][QIN.numberOfBCnodes] = ON[BSE ];
                  Q.q27[BNW ][QIN.numberOfBCnodes] = ON[BNW ];

                  QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void findKforQ(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

   VF_LOG_CRITICAL("findKforQ() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m, mm, l;
   real test = (real)0.f;
   int nx                       = para->getParH(lev)->nx;
   int ny                       = para->getParH(lev)->ny;
   unsigned int nnx             = para->getParH(lev)->gridNX;
   unsigned int nny             = para->getParH(lev)->gridNY;
   unsigned int nnz             = para->getParH(lev)->gridNZ;
   int* geo_mat                 = para->getParH(lev)->geo;
   QforBoundaryConditions &QIN  = para->getParH(lev)->noSlipBC;
   QIN.numberOfBCnodes = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   for(k=STARTOFFZ + 1 ; k<=nnz+STARTOFFZ-2 ; k++){
      for(j=STARTOFFY + 1 ; j<=nny+STARTOFFY-2 ; j++){       //j<=nny/2+STARTOFFY   //j<=STARTOFFY+1
         for(i=STARTOFFX + 1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test =(real)0.f;
               for(l=0;l<=26;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if((geo_mat[mm] == GEO_SOLID) || (geo_mat[mm] == GEO_VOID)){
                     ON[l] =(real) 1.f;
                  }
                  else{
                     ON[l] = (real)0.f;
                  }
                  test += ON[l];
               }
               if (test>0)
               {
                  QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void findQ_MG( int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk, unsigned int sizeQ, real* QQ, QforBoundaryConditions &QIN)
{
   QforBoundaryConditions Q;
   Q.q27[E   ] = &QQ[E   *sizeQ];
   Q.q27[W   ] = &QQ[W   *sizeQ];
   Q.q27[N   ] = &QQ[N   *sizeQ];
   Q.q27[S   ] = &QQ[S   *sizeQ];
   Q.q27[T   ] = &QQ[T   *sizeQ];
   Q.q27[B   ] = &QQ[B   *sizeQ];
   Q.q27[NE  ] = &QQ[NE  *sizeQ];
   Q.q27[SW  ] = &QQ[SW  *sizeQ];
   Q.q27[SE  ] = &QQ[SE  *sizeQ];
   Q.q27[NW  ] = &QQ[NW  *sizeQ];
   Q.q27[TE  ] = &QQ[TE  *sizeQ];
   Q.q27[BW  ] = &QQ[BW  *sizeQ];
   Q.q27[BE  ] = &QQ[BE  *sizeQ];
   Q.q27[TW  ] = &QQ[TW  *sizeQ];
   Q.q27[TN  ] = &QQ[TN  *sizeQ];
   Q.q27[BS  ] = &QQ[BS  *sizeQ];
   Q.q27[BN  ] = &QQ[BN  *sizeQ];
   Q.q27[TS  ] = &QQ[TS  *sizeQ];
   Q.q27[REST] = &QQ[REST*sizeQ];
   Q.q27[TNE ] = &QQ[TNE *sizeQ];
   Q.q27[TSW ] = &QQ[TSW *sizeQ];
   Q.q27[TSE ] = &QQ[TSE *sizeQ];
   Q.q27[TNW ] = &QQ[TNW *sizeQ];
   Q.q27[BNE ] = &QQ[BNE *sizeQ];
   Q.q27[BSW ] = &QQ[BSW *sizeQ];
   Q.q27[BSE ] = &QQ[BSE *sizeQ];
   Q.q27[BNW ] = &QQ[BNW *sizeQ];

   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

    VF_LOG_CRITICAL("findQ_MG() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   unsigned int i, j, k, m, mm, l;
   int relx, rely, relz;

   unsigned int centerX = nnx / 2;
   unsigned int centerY = nny / 2;
   unsigned int centerZ = nnz / 2;
   real      radius  = nny / 2.56f;

   QIN.numberOfBCnodes = 0;

   for(k=STARTOFFZ+2 ; k<=nnz+STARTOFFZ-3 ; k++){
      for(j=STARTOFFY+2 ; j<=nny+STARTOFFY-3 ; j++){
         for(i=STARTOFFX +2; i<=nnx+STARTOFFX-3 ; i++){
            m = nx*(ny*k + j) + i;
            if((geo_mat[m] == GEO_SOLID) || (geo_mat[m] == GEO_VOID)){
               relx = i - STARTOFFX - centerX;
               rely = j - STARTOFFY - centerY;
               relz = k - STARTOFFZ - centerZ;
               ON[18] = (real)-1.f;
               for(l=0;l<=26;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if(geo_mat[mm] == GEO_FLUID){
                     ON[l] = -(((real)ex[l]*(real)relx + (real)ey[l]*(real)rely + (real)ez[l]*(real)relz +
                        sqrt(pow((real)ex[l]*(real)relx + (real)ey[l]*(real)rely + (real)ez[l]*(real)relz,2) +
                        (pow((real)ex[l],2) + pow((real)ey[l],2) + pow((real)ez[l],2))*
                        (pow(radius,2) - pow((real)relx,2) - pow((real)rely,2) -
                        pow((real)relz,2))))/(pow((real)ex[l],2) + pow((real)ey[l],2) + pow((real)ez[l],2)));

                     ON[18] = (real)1.f; //ZERO
                  }
                  else{
                     ON[l] = (real)-1.f;
                  }
               }
               if (ON[18]==1.f)
               {
                  QIN.k[QIN.numberOfBCnodes]          = kk[m];

                  Q.q27[E   ][QIN.numberOfBCnodes] = ON[E   ];
                  Q.q27[W   ][QIN.numberOfBCnodes] = ON[W   ];
                  Q.q27[N   ][QIN.numberOfBCnodes] = ON[N   ];
                  Q.q27[S   ][QIN.numberOfBCnodes] = ON[S   ];
                  Q.q27[T   ][QIN.numberOfBCnodes] = ON[T   ];
                  Q.q27[B   ][QIN.numberOfBCnodes] = ON[B   ];
                  Q.q27[NE  ][QIN.numberOfBCnodes] = ON[NE  ];
                  Q.q27[SW  ][QIN.numberOfBCnodes] = ON[SW  ];
                  Q.q27[SE  ][QIN.numberOfBCnodes] = ON[SE  ];
                  Q.q27[NW  ][QIN.numberOfBCnodes] = ON[NW  ];
                  Q.q27[TE  ][QIN.numberOfBCnodes] = ON[TE  ];
                  Q.q27[BW  ][QIN.numberOfBCnodes] = ON[BW  ];
                  Q.q27[BE  ][QIN.numberOfBCnodes] = ON[BE  ];
                  Q.q27[TW  ][QIN.numberOfBCnodes] = ON[TW  ];
                  Q.q27[TN  ][QIN.numberOfBCnodes] = ON[TN  ];
                  Q.q27[BS  ][QIN.numberOfBCnodes] = ON[BS  ];
                  Q.q27[BN  ][QIN.numberOfBCnodes] = ON[BN  ];
                  Q.q27[TS  ][QIN.numberOfBCnodes] = ON[TS  ];
                  Q.q27[REST][QIN.numberOfBCnodes] = ON[REST];
                  Q.q27[TNE ][QIN.numberOfBCnodes] = ON[TNE ];
                  Q.q27[TSW ][QIN.numberOfBCnodes] = ON[TSW ];
                  Q.q27[TSE ][QIN.numberOfBCnodes] = ON[TSE ];
                  Q.q27[TNW ][QIN.numberOfBCnodes] = ON[TNW ];
                  Q.q27[BNE ][QIN.numberOfBCnodes] = ON[BNE ];
                  Q.q27[BSW ][QIN.numberOfBCnodes] = ON[BSW ];
                  Q.q27[BSE ][QIN.numberOfBCnodes] = ON[BSE ];
                  Q.q27[BNW ][QIN.numberOfBCnodes] = ON[BNW ];

                  QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void findKforQ_MG(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, QforBoundaryConditions &QIN)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

    VF_LOG_CRITICAL("findKforQ_MG() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   unsigned int i, j, k, m, mm, l;
   real test = (real)0.f;

   QIN.numberOfBCnodes = 0;

   for(k=STARTOFFZ+2 ; k<=nnz+STARTOFFZ-3 ; k++){
      for(j=STARTOFFY+2 ; j<=nny+STARTOFFY-3 ; j++){
         for(i=STARTOFFX +2; i<=nnx+STARTOFFX-3 ; i++){
            m = nx*(ny*k + j) + i;
            if((geo_mat[m] == GEO_SOLID) || (geo_mat[m] == GEO_VOID)){
               test =(real)0.f;
               for(l=0;l<=26;l++){
                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if(geo_mat[mm] == GEO_FLUID){
                     ON[l] = (real)1.f;
                  }
                  else{
                     ON[l] = (real)0.f;
                  }
                  test += ON[l];
               }
               if (test>0)
               {
                  QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
void findQInflow(Parameter* para)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findQInflow() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m;//, mm, l;
   //int nx                        = para->getParH(para->getFine())->nx;
   //int ny                        = para->getParH(para->getFine())->ny;
   //unsigned int nnx              = para->getParH(para->getFine())->gridNX;
   //unsigned int nny              = para->getParH(para->getFine())->gridNY;
   //unsigned int nnz              = para->getParH(para->getFine())->gridNZ;
   //int* geo_mat                  = para->getParH(para->getFine())->geo;
   //unsigned int* kk              = para->getParH(para->getFine())->k;
   int nx                        = para->getParH(para->getCoarse())->nx;
   int ny                        = para->getParH(para->getCoarse())->ny;
   unsigned int nnx              = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny              = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz              = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                  = para->getParH(para->getCoarse())->geo;
   unsigned int* kk              = para->getParH(para->getCoarse())->k;
   unsigned int sizeQ            = para->getParH(para->getCoarse())->velocityBC.numberOfBCnodes;
   //real* rhoBC                = para->getParH(para->getCoarse())->velocityBC.RhoBC;
   real u0                    = para->getVelocity();
   real* vx                   = para->getParH(para->getCoarse())->velocityBC.Vx;
   real* vy                   = para->getParH(para->getCoarse())->velocityBC.Vy;
   real* vz                   = para->getParH(para->getCoarse())->velocityBC.Vz;
   real*deltaVz               = para->getParH(para->getCoarse())->velocityBC.deltaVz;
   real* QQ                   = para->getParH(para->getCoarse())->velocityBC.q27[0];
   QforBoundaryConditions &QIN   = para->getParH(para->getCoarse())->velocityBC;
   //unsigned int nxny = nx*ny;
   QIN.numberOfBCnodes = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   QforBoundaryConditions Q;
   Q.q27[E   ] = &QQ[E   *sizeQ];
   Q.q27[W   ] = &QQ[W   *sizeQ];
   Q.q27[N   ] = &QQ[N   *sizeQ];
   Q.q27[S   ] = &QQ[S   *sizeQ];
   Q.q27[T   ] = &QQ[T   *sizeQ];
   Q.q27[B   ] = &QQ[B   *sizeQ];
   Q.q27[NE  ] = &QQ[NE  *sizeQ];
   Q.q27[SW  ] = &QQ[SW  *sizeQ];
   Q.q27[SE  ] = &QQ[SE  *sizeQ];
   Q.q27[NW  ] = &QQ[NW  *sizeQ];
   Q.q27[TE  ] = &QQ[TE  *sizeQ];
   Q.q27[BW  ] = &QQ[BW  *sizeQ];
   Q.q27[BE  ] = &QQ[BE  *sizeQ];
   Q.q27[TW  ] = &QQ[TW  *sizeQ];
   Q.q27[TN  ] = &QQ[TN  *sizeQ];
   Q.q27[BS  ] = &QQ[BS  *sizeQ];
   Q.q27[BN  ] = &QQ[BN  *sizeQ];
   Q.q27[TS  ] = &QQ[TS  *sizeQ];
   Q.q27[REST] = &QQ[REST*sizeQ];
   Q.q27[TNE ] = &QQ[TNE *sizeQ];
   Q.q27[TSW ] = &QQ[TSW *sizeQ];
   Q.q27[TSE ] = &QQ[TSE *sizeQ];
   Q.q27[TNW ] = &QQ[TNW *sizeQ];
   Q.q27[BNE ] = &QQ[BNE *sizeQ];
   Q.q27[BSW ] = &QQ[BSW *sizeQ];
   Q.q27[BSE ] = &QQ[BSE *sizeQ];
   Q.q27[BNW ] = &QQ[BNW *sizeQ];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //unsigned int li = ((nnx+STARTOFFX-2)-(STARTOFFX+1)-1);
   //unsigned int lj = ((nny+STARTOFFY-2)-(STARTOFFY+1)-1);

   //k = STARTOFFZ + 1;
   k = nnz+STARTOFFZ - 1/*3*/;

   for(j=STARTOFFY/*+1*/ ; j<=nny+STARTOFFY/*-2*/ ; j++){
       for(i=STARTOFFX/*+1*/; i<=nnx+STARTOFFX/*-2*/ ; i++){
         m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               QIN.k[QIN.numberOfBCnodes]          = kk[m];
               //vx[QIN.numberOfBCnodes]             = (real)0.f;
			   vx[QIN.numberOfBCnodes]             = u0;
               vy[QIN.numberOfBCnodes]             = (real)0.f;
               vz[QIN.numberOfBCnodes]             = (real)0.f;
               //vz[QIN.numberOfBCnodes]             = u0;
               //vz[QIN.numberOfBCnodes]             = (real)(u0*2.f)*((-4.f*i*i + nnx*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*nnx + 8.f*STARTOFFX))*(-4.f*j*j + nny*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*nny + 8.f*STARTOFFY)))/((2.f - nnx)*(2.f - nnx)*(2.f - nny)*(2.f - nny));
               //vz[QIN.numberOfBCnodes]             = (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);
               //vz[QIN.numberOfBCnodes]             = (real)(16.f*(u0*2.f)*i*j*(nx-i)*(ny-j))/(nx*nx*ny*ny);
               deltaVz[QIN.numberOfBCnodes]        = (real)0.f;
			   //////////////////////////////////////////////////////////////////////////
               //Q.q27[E   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[W   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[N   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[S   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[T   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[B   ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[NE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[SW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[SE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[NW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[BW  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[BE  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[TW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TN  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[BS  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[BN  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[TS  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[REST][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TNE ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TSW ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TSE ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[TNW ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[BNE ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[BSW ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[BSE ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[BNW ][QIN.numberOfBCnodes] = (real)1.f;
			   //////////////////////////////////////////////////////////////////////////


               // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

			   Q.q27[E   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[W   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[N   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[S   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[T   ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[B   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[NE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[SW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[SE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[NW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[TE  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[BW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[BE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[TW  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[TN  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[BS  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[BN  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[TS  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[REST][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[TNE ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[TSW ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[TSE ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[TNW ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[BNE ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[BSW ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[BSE ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[BNW ][QIN.numberOfBCnodes] = (real)-1.f;
			   //////////////////////////////////////////////////////////////////////////
			   QIN.numberOfBCnodes++;
            }
       }
   }
   //k = STARTOFFZ +2;
   //{
   //   for(j=STARTOFFY ; j<=nny+STARTOFFY-1 ; j++){
   //      for(i=STARTOFFX; i<=nnx+STARTOFFX-1 ; i++){
   //         m = nx*(ny*k + j) + i;
   //         if(geo_mat[m]==GEO_FLUID){
   //            ON[18] = -1.f;
   //            for(l=0;l<=26;l++){
   //               mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
   //               if(ez[l]==-1){
   //                  ON[l] = 1.f;

   //                  ON[18] = 1.f; //ZERO
   //               }
   //               else{
   //                  ON[l] = -1.f;
   //               }
   //            }
   //            if (ON[18]==1.f)
   //            {
   //               QIN.k[QIN.numberOfBCnodes]          = m;
   //               vx[QIN.numberOfBCnodes]             = 0.f;
   //               vy[QIN.numberOfBCnodes]             = 0.f;
   //               vz[QIN.numberOfBCnodes]             = u0;

   //               Q.q27[E   ][QIN.numberOfBCnodes] = ON[E   ];
   //               Q.q27[W   ][QIN.numberOfBCnodes] = ON[W   ];
   //               Q.q27[N   ][QIN.numberOfBCnodes] = ON[N   ];
   //               Q.q27[S   ][QIN.numberOfBCnodes] = ON[S   ];
   //               Q.q27[T   ][QIN.numberOfBCnodes] = ON[T   ];
   //               Q.q27[B   ][QIN.numberOfBCnodes] = ON[B   ];
   //               Q.q27[NE  ][QIN.numberOfBCnodes] = ON[NE  ];
   //               Q.q27[SW  ][QIN.numberOfBCnodes] = ON[SW  ];
   //               Q.q27[SE  ][QIN.numberOfBCnodes] = ON[SE  ];
   //               Q.q27[NW  ][QIN.numberOfBCnodes] = ON[NW  ];
   //               Q.q27[TE  ][QIN.numberOfBCnodes] = ON[TE  ];
   //               Q.q27[BW  ][QIN.numberOfBCnodes] = ON[BW  ];
   //               Q.q27[BE  ][QIN.numberOfBCnodes] = ON[BE  ];
   //               Q.q27[TW  ][QIN.numberOfBCnodes] = ON[TW  ];
   //               Q.q27[TN  ][QIN.numberOfBCnodes] = ON[TN  ];
   //               Q.q27[BS  ][QIN.numberOfBCnodes] = ON[BS  ];
   //               Q.q27[BN  ][QIN.numberOfBCnodes] = ON[BN  ];
   //               Q.q27[TS  ][QIN.numberOfBCnodes] = ON[TS  ];
   //               Q.q27[REST][QIN.numberOfBCnodes] = ON[REST];
   //               Q.q27[TNE ][QIN.numberOfBCnodes] = ON[TNE ];
   //               Q.q27[TSW ][QIN.numberOfBCnodes] = ON[TSW ];
   //               Q.q27[TSE ][QIN.numberOfBCnodes] = ON[TSE ];
   //               Q.q27[TNW ][QIN.numberOfBCnodes] = ON[TNW ];
   //               Q.q27[BNE ][QIN.numberOfBCnodes] = ON[BNE ];
   //               Q.q27[BSW ][QIN.numberOfBCnodes] = ON[BSW ];
   //               Q.q27[BSE ][QIN.numberOfBCnodes] = ON[BSE ];
   //               Q.q27[BNW ][QIN.numberOfBCnodes] = ON[BNW ];

   //               QIN.numberOfBCnodes++;
   //            }
   //         }
   //      }
   //   }
   //}

   //for(k=STARTOFFZ ; k<=nnz+STARTOFFZ-1 ; k++){
   //   for(j=STARTOFFY ; j<=nny+STARTOFFY-1 ; j++){
   //      for(i=STARTOFFX; i<=nnx+STARTOFFX-1 ; i++){
   //         m = nx*(ny*k + j) + i;
   //         if(geo_mat[m]==GEO_FLUID){
   //            ON[18] = -1.f;
   //            for(l=0;l<=26;l++){
   //               mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
   //               if(geo_mat[mm] == GEO_VELO){
   //                  ON[l] = 1.f;

   //                  ON[18] = 1.f; //ZERO
   //               }
   //               else{
   //                  ON[l] = -1.f;
   //               }
   //            }
   //            if (ON[18]==1.f)
   //            {
   //               QIN.k[QIN.numberOfBCnodes]          = m;
   //               vx[QIN.numberOfBCnodes]             = u0;//0.f;
   //               vy[QIN.numberOfBCnodes]             = 0.f;
   //               vz[QIN.numberOfBCnodes]             = 0.f;//u0;

   //               Q.q27[E   ][QIN.numberOfBCnodes] = ON[E   ];
   //               Q.q27[W   ][QIN.numberOfBCnodes] = ON[W   ];
   //               Q.q27[N   ][QIN.numberOfBCnodes] = ON[N   ];
   //               Q.q27[S   ][QIN.numberOfBCnodes] = ON[S   ];
   //               Q.q27[T   ][QIN.numberOfBCnodes] = ON[T   ];
   //               Q.q27[B   ][QIN.numberOfBCnodes] = ON[B   ];
   //               Q.q27[NE  ][QIN.numberOfBCnodes] = ON[NE  ];
   //               Q.q27[SW  ][QIN.numberOfBCnodes] = ON[SW  ];
   //               Q.q27[SE  ][QIN.numberOfBCnodes] = ON[SE  ];
   //               Q.q27[NW  ][QIN.numberOfBCnodes] = ON[NW  ];
   //               Q.q27[TE  ][QIN.numberOfBCnodes] = ON[TE  ];
   //               Q.q27[BW  ][QIN.numberOfBCnodes] = ON[BW  ];
   //               Q.q27[BE  ][QIN.numberOfBCnodes] = ON[BE  ];
   //               Q.q27[TW  ][QIN.numberOfBCnodes] = ON[TW  ];
   //               Q.q27[TN  ][QIN.numberOfBCnodes] = ON[TN  ];
   //               Q.q27[BS  ][QIN.numberOfBCnodes] = ON[BS  ];
   //               Q.q27[BN  ][QIN.numberOfBCnodes] = ON[BN  ];
   //               Q.q27[TS  ][QIN.numberOfBCnodes] = ON[TS  ];
   //               Q.q27[REST][QIN.numberOfBCnodes] = ON[REST];
   //               Q.q27[TNE ][QIN.numberOfBCnodes] = ON[TNE ];
   //               Q.q27[TSW ][QIN.numberOfBCnodes] = ON[TSW ];
   //               Q.q27[TSE ][QIN.numberOfBCnodes] = ON[TSE ];
   //               Q.q27[TNW ][QIN.numberOfBCnodes] = ON[TNW ];
   //               Q.q27[BNE ][QIN.numberOfBCnodes] = ON[BNE ];
   //               Q.q27[BSW ][QIN.numberOfBCnodes] = ON[BSW ];
   //               Q.q27[BSE ][QIN.numberOfBCnodes] = ON[BSE ];
   //               Q.q27[BNW ][QIN.numberOfBCnodes] = ON[BNW ];

   //               QIN.numberOfBCnodes++;
   //            }
   //         }
   //      }
   //   }
   //}
}

////////////////////////////////////////////////////////////////////////////////
void findKforQInflow(Parameter* para)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQInflow() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   //int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   //int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //unsigned int mm;
   unsigned int i, j, k, m, l;
   real test = 0.f;
   //int nx                        = para->getParH(para->getFine())->nx;
   //int ny                        = para->getParH(para->getFine())->ny;
   //unsigned int nnx              = para->getParH(para->getFine())->gridNX;
   //unsigned int nny              = para->getParH(para->getFine())->gridNY;
   //unsigned int nnz              = para->getParH(para->getFine())->gridNZ;
   //int* geo_mat                  = para->getParH(para->getFine())->geo;
   int nx                        = para->getParH(para->getCoarse())->nx;
   int ny                        = para->getParH(para->getCoarse())->ny;
   unsigned int nnx              = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny              = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz              = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                  = para->getParH(para->getCoarse())->geo;
   QforBoundaryConditions &QIN   = para->getParH(para->getCoarse())->velocityBC;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   QIN.numberOfBCnodes = 0;

   k=nnz+STARTOFFZ-1/*3*/;
   //k=STARTOFFZ+1;
   {
      for(j=STARTOFFY/*+1*/ ; j<=nny+STARTOFFY/*-2*/ ; j++){
         for(i=STARTOFFX/*+1*/; i<=nnx+STARTOFFX/*-2*/ ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test = (real)0.f;
               for(l=0;l<=26;l++){
                  //mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if(ez[l]==1/*-1*/){
                     ON[l] = (real) 1.f;
                  }
                  else{
                     ON[l] = (real) 0.f;
                  }
                  test += ON[l];
               }
               if (test>0)
               {
				   QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }

   //for(k=STARTOFFZ ; k<=nnz+STARTOFFZ-1 ; k++){
   //   for(j=STARTOFFY ; j<=nny+STARTOFFY-1 ; j++){
   //      for(i=STARTOFFX; i<=nnx+STARTOFFX-1 ; i++){
   //         m = nx*(ny*k + j) + i;
   //         if(geo_mat[m]==GEO_FLUID){
   //            test =0.f;
   //            for(l=0;l<=26;l++){
   //               mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
   //               if(geo_mat[mm] == GEO_VELO){
   //                  ON[l] = 1.f;
   //               }
   //               else{
   //                  ON[l] = 0.f;
   //               }
   //               test += ON[l];
   //            }
   //            if (test>0)
   //            {
   //               QIN.numberOfBCnodes++;
   //            }
   //         }
   //      }
   //   }
   //}
}


////////////////////////////////////////////////////////////////////////////////
void findQOutflow(Parameter* para)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findQOutflow() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   //int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   //int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   //int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   //real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   unsigned int i, j, k, m;//, mm, l;
   //int nx                        = para->getParH(para->getFine())->nx;
   //int ny                        = para->getParH(para->getFine())->ny;
   //unsigned int nnx              = para->getParH(para->getFine())->gridNX;
   //unsigned int nny              = para->getParH(para->getFine())->gridNY;
   //unsigned int nnz              = para->getParH(para->getFine())->gridNZ;
   //int* geo_mat                  = para->getParH(para->getFine())->geo;
   //unsigned int* kk              = para->getParH(para->getFine())->k;
   int nx                        = para->getParH(para->getCoarse())->nx;
   int ny                        = para->getParH(para->getCoarse())->ny;
   unsigned int nnx              = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny              = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz              = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                  = para->getParH(para->getCoarse())->geo;
   unsigned int* kk              = para->getParH(para->getCoarse())->k;
   unsigned int sizeQ            = para->getParH(para->getCoarse())->outflowBC.numberOfBCnodes;
   real* rhoBC                = para->getParH(para->getCoarse())->outflowBC.RhoBC;
   real u0                    = para->getVelocity();
   real* vx                   = para->getParH(para->getCoarse())->outflowBC.Vx;
   real* vy                   = para->getParH(para->getCoarse())->outflowBC.Vy;
   real* vz                   = para->getParH(para->getCoarse())->outflowBC.Vz;
   real*deltaVz               = para->getParH(para->getCoarse())->outflowBC.deltaVz;
   real* QQ                   = para->getParH(para->getCoarse())->outflowBC.q27[0];
   QforBoundaryConditions &QIN   = para->getParH(para->getCoarse())->outflowBC;
   unsigned int nxny = nx*ny;
   QIN.numberOfBCnodes = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   QforBoundaryConditions Q;
   Q.q27[E   ] = &QQ[E   *sizeQ];
   Q.q27[W   ] = &QQ[W   *sizeQ];
   Q.q27[N   ] = &QQ[N   *sizeQ];
   Q.q27[S   ] = &QQ[S   *sizeQ];
   Q.q27[T   ] = &QQ[T   *sizeQ];
   Q.q27[B   ] = &QQ[B   *sizeQ];
   Q.q27[NE  ] = &QQ[NE  *sizeQ];
   Q.q27[SW  ] = &QQ[SW  *sizeQ];
   Q.q27[SE  ] = &QQ[SE  *sizeQ];
   Q.q27[NW  ] = &QQ[NW  *sizeQ];
   Q.q27[TE  ] = &QQ[TE  *sizeQ];
   Q.q27[BW  ] = &QQ[BW  *sizeQ];
   Q.q27[BE  ] = &QQ[BE  *sizeQ];
   Q.q27[TW  ] = &QQ[TW  *sizeQ];
   Q.q27[TN  ] = &QQ[TN  *sizeQ];
   Q.q27[BS  ] = &QQ[BS  *sizeQ];
   Q.q27[BN  ] = &QQ[BN  *sizeQ];
   Q.q27[TS  ] = &QQ[TS  *sizeQ];
   Q.q27[REST] = &QQ[REST*sizeQ];
   Q.q27[TNE ] = &QQ[TNE *sizeQ];
   Q.q27[TSW ] = &QQ[TSW *sizeQ];
   Q.q27[TSE ] = &QQ[TSE *sizeQ];
   Q.q27[TNW ] = &QQ[TNW *sizeQ];
   Q.q27[BNE ] = &QQ[BNE *sizeQ];
   Q.q27[BSW ] = &QQ[BSW *sizeQ];
   Q.q27[BSE ] = &QQ[BSE *sizeQ];
   Q.q27[BNW ] = &QQ[BNW *sizeQ];


   //unsigned int li = ((nnx+STARTOFFX-2)-(STARTOFFX+1)-1);
   //unsigned int lj = ((nny+STARTOFFY-2)-(STARTOFFY+1)-1);

   k = nnz + STARTOFFZ - 3;

   for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
       for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
         m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               QIN.k[QIN.numberOfBCnodes]          = kk[m];
               QIN.kN[QIN.numberOfBCnodes]         = kk[m-nxny];
               rhoBC[QIN.numberOfBCnodes]          = (real)0.f;
               vx[QIN.numberOfBCnodes]             = (real)0.f;
               vy[QIN.numberOfBCnodes]             = (real)0.f;
			   //vz[QIN.numberOfBCnodes]             = u0;
               vz[QIN.numberOfBCnodes]             = (real)(u0*2.f)*((-4.f*i*i + nnx*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*nnx + 8.f*STARTOFFX))*(-4.f*j*j + nny*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*nny + 8.f*STARTOFFY)))/((2.f - nnx)*(2.f - nnx)*(2.f - nny)*(2.f - nny));
               //vz[QIN.numberOfBCnodes]             =  (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);
               //vz[QIN.numberOfBCnodes]             = (real)(16.f*(u0*2.f)*i*j*(nx-i)*(ny-j))/(nx*nx*ny*ny);
               deltaVz[QIN.numberOfBCnodes]        = (real)0.f;
               Q.q27[E   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[W   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[N   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[S   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[T   ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[B   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[NE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[SW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[SE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[NW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[TE  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[BW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[BE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[TW  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[TN  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[BS  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[BN  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[TS  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[REST][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[TNE ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[TSW ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[TSE ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[TNW ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[BNE ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[BSW ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[BSE ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[BNW ][QIN.numberOfBCnodes] = (real)-1.f;
               QIN.numberOfBCnodes++;
            }
       }
   }

   //i = nnx / 2 + STARTOFFX;
   //j = nny / 2 + STARTOFFY;
   //k = nnz / 2 + STARTOFFZ;
   //QIN.numberOfBCnodes = 0;
   //rhoBC[QIN.numberOfBCnodes]        = 0.1f;
   //QIN.numberOfBCnodes++;

}

////////////////////////////////////////////////////////////////////////////////
void findKforQOutflow(Parameter* para)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQOutflow() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
   //int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
   //int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
   real ON[27];
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //unsigned int mm;
   unsigned int i, j, k, m, l;
   real test = (real) 0.f;
   //int nx                        = para->getParH(para->getFine())->nx;
   //int ny                        = para->getParH(para->getFine())->ny;
   //unsigned int nnx              = para->getParH(para->getFine())->gridNX;
   //unsigned int nny              = para->getParH(para->getFine())->gridNY;
   //unsigned int nnz              = para->getParH(para->getFine())->gridNZ;
   //int* geo_mat                  = para->getParH(para->getFine())->geo;
   int nx                        = para->getParH(para->getCoarse())->nx;
   int ny                        = para->getParH(para->getCoarse())->ny;
   unsigned int nnx              = para->getParH(para->getCoarse())->gridNX;
   unsigned int nny              = para->getParH(para->getCoarse())->gridNY;
   unsigned int nnz              = para->getParH(para->getCoarse())->gridNZ;
   int* geo_mat                  = para->getParH(para->getCoarse())->geo;
   QforBoundaryConditions &QIN   = para->getParH(para->getCoarse())->outflowBC;
   QIN.numberOfBCnodes = 0;
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   k=nnz+STARTOFFZ-3;
   {
      for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
            m = nx*(ny*k + j) + i;
            if(geo_mat[m]==GEO_FLUID){
               test = (real)0.f;
               for(l=0;l<=26;l++){
                  //mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
                  if(ez[l]==1){
                     ON[l] = (real) 1.f;
                  }
                  else{
                     ON[l] = (real) 0.f;
                  }
                  test += ON[l];
               }
               if (test>0)
               {
                  QIN.numberOfBCnodes++;
               }
            }
         }
      }
   }

   //QIN.numberOfBCnodes = 0;
   //QIN.numberOfBCnodes++;

}
// TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
//////////////////////////////////////////////////////////////////////////////////
//void findQSchlaff( int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat, unsigned int* kk,
//                   unsigned int sizeQN, real* vxN, real* vyN, real* vzN, real*deltaVN, real* QQN, QforBoundaryConditions &QNin,
//                   unsigned int sizeQS, real* vxS, real* vyS, real* vzS, real*deltaVS, real* QQS, QforBoundaryConditions &QSin,
//                   unsigned int sizeQE, real* vxE, real* vyE, real* vzE, real*deltaVE, real* QQE, QforBoundaryConditions &QEin,
//                   unsigned int sizeQW, real* vxW, real* vyW, real* vzW, real*deltaVW, real* QQW, QforBoundaryConditions &QWin)
//{
//   QforBoundaryConditions QN;
//   QN.q27[E   ] = &QQN[E   *sizeQN];
//   QN.q27[W   ] = &QQN[W   *sizeQN];
//   QN.q27[N   ] = &QQN[N   *sizeQN];
//   QN.q27[S   ] = &QQN[S   *sizeQN];
//   QN.q27[T   ] = &QQN[T   *sizeQN];
//   QN.q27[B   ] = &QQN[B   *sizeQN];
//   QN.q27[NE  ] = &QQN[NE  *sizeQN];
//   QN.q27[SW  ] = &QQN[SW  *sizeQN];
//   QN.q27[SE  ] = &QQN[SE  *sizeQN];
//   QN.q27[NW  ] = &QQN[NW  *sizeQN];
//   QN.q27[TE  ] = &QQN[TE  *sizeQN];
//   QN.q27[BW  ] = &QQN[BW  *sizeQN];
//   QN.q27[BE  ] = &QQN[BE  *sizeQN];
//   QN.q27[TW  ] = &QQN[TW  *sizeQN];
//   QN.q27[TN  ] = &QQN[TN  *sizeQN];
//   QN.q27[BS  ] = &QQN[BS  *sizeQN];
//   QN.q27[BN  ] = &QQN[BN  *sizeQN];
//   QN.q27[TS  ] = &QQN[TS  *sizeQN];
//   QN.q27[REST] = &QQN[REST*sizeQN];
//   QN.q27[TNE ] = &QQN[TNE *sizeQN];
//   QN.q27[TSW ] = &QQN[TSW *sizeQN];
//   QN.q27[TSE ] = &QQN[TSE *sizeQN];
//   QN.q27[TNW ] = &QQN[TNW *sizeQN];
//   QN.q27[BNE ] = &QQN[BNE *sizeQN];
//   QN.q27[BSW ] = &QQN[BSW *sizeQN];
//   QN.q27[BSE ] = &QQN[BSE *sizeQN];
//   QN.q27[BNW ] = &QQN[BNW *sizeQN];
//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
//   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
//   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
//   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
//   //real ON[27];
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//}
//
//////////////////////////////////////////////////////////////////////////////////
//void findKforQSchlaff(int nx, int ny, unsigned int nnx, unsigned int nny, unsigned int nnz, int* geo_mat,
//                      QforBoundaryConditions &QN, QforBoundaryConditions &QS, QforBoundaryConditions &QE, QforBoundaryConditions &QW)
//{
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
//   int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
//   int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
//   int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
//   real ON[27];
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   unsigned int i, j, k, m, mm, l;
//   real test = (real) 0.f;
//
//   QN.numberOfBCnodes = 0;
//   QS.numberOfBCnodes = 0;
//   QE.numberOfBCnodes = 0;
//   QW.numberOfBCnodes = 0;
//
//   for(k=STARTOFFZ+1 ; k<=nnz+STARTOFFZ-2 ; k++){
//      for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
//         for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
//            m = nx*(ny*k + j) + i;
//            if(geo_mat[m]==GEO_FLUID){
//               test = (real)0.f;
//               for(l=0;l<=26;l++){
//                  mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
//                  if(ez[l]==1){
//                     ON[l] = (real) 1.f;
//                  }
//                  else{
//                     ON[l] = (real) 0.f;
//                  }
//                  test += ON[l];
//               }
//               if (test>0)
//               {
//                  QN.numberOfBCnodes++;
//               }
//            }
//         }
//      }
//   }
//}



////////////////////////////////////////////////////////////////////////////////
void findQPressX0(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQPressX0() is deprecated! - see comment above for more information");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
	//int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
	//int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
	//int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int i, j, k, m;
	int nx                        = para->getParH(lev)->nx;
	int ny                        = para->getParH(lev)->ny;
	unsigned int nnx              = para->getParH(lev)->gridNX;
	unsigned int nny              = para->getParH(lev)->gridNY;
	unsigned int nnz              = para->getParH(lev)->gridNZ;
	int* geo_mat                  = para->getParH(lev)->geo;
	unsigned int* kk              = para->getParH(lev)->k;
	//unsigned int sizeQ            = para->getParH(lev)->outflowBC.numberOfBCnodes;
	unsigned int sizeQ            = para->getParH(lev)->QpressX0.numberOfBCnodes;
	real* rhoBC                = para->getParH(lev)->QpressX0.RhoBC;
	real u0                    = para->getVelocity();
	real* vx                   = para->getParH(lev)->QpressX0.Vx;
	real* vy                   = para->getParH(lev)->QpressX0.Vy;
	real* vz                   = para->getParH(lev)->QpressX0.Vz;
	real*deltaVz               = para->getParH(lev)->QpressX0.deltaVz;
	real* QQ                   = para->getParH(lev)->QpressX0.q27[0];
	QforBoundaryConditions &QIN   = para->getParH(lev)->QpressX0;
	//unsigned int nxny = nx*ny;
	QIN.numberOfBCnodes = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QforBoundaryConditions Q;
	Q.q27[E   ] = &QQ[E   *sizeQ];
	Q.q27[W   ] = &QQ[W   *sizeQ];
	Q.q27[N   ] = &QQ[N   *sizeQ];
	Q.q27[S   ] = &QQ[S   *sizeQ];
	Q.q27[T   ] = &QQ[T   *sizeQ];
	Q.q27[B   ] = &QQ[B   *sizeQ];
	Q.q27[NE  ] = &QQ[NE  *sizeQ];
	Q.q27[SW  ] = &QQ[SW  *sizeQ];
	Q.q27[SE  ] = &QQ[SE  *sizeQ];
	Q.q27[NW  ] = &QQ[NW  *sizeQ];
	Q.q27[TE  ] = &QQ[TE  *sizeQ];
	Q.q27[BW  ] = &QQ[BW  *sizeQ];
	Q.q27[BE  ] = &QQ[BE  *sizeQ];
	Q.q27[TW  ] = &QQ[TW  *sizeQ];
	Q.q27[TN  ] = &QQ[TN  *sizeQ];
	Q.q27[BS  ] = &QQ[BS  *sizeQ];
	Q.q27[BN  ] = &QQ[BN  *sizeQ];
	Q.q27[TS  ] = &QQ[TS  *sizeQ];
	Q.q27[REST] = &QQ[REST*sizeQ];
	Q.q27[TNE ] = &QQ[TNE *sizeQ];
	Q.q27[TSW ] = &QQ[TSW *sizeQ];
	Q.q27[TSE ] = &QQ[TSE *sizeQ];
	Q.q27[TNW ] = &QQ[TNW *sizeQ];
	Q.q27[BNE ] = &QQ[BNE *sizeQ];
	Q.q27[BSW ] = &QQ[BSW *sizeQ];
	Q.q27[BSE ] = &QQ[BSE *sizeQ];
	Q.q27[BNW ] = &QQ[BNW *sizeQ];


	//unsigned int li = ((nnx+STARTOFFX-2)-(STARTOFFX+1)-1);
	//unsigned int lj = ((nny+STARTOFFY-2)-(STARTOFFY+1)-1);

	i=STARTOFFX+1;
	//k=nnz+STARTOFFZ-3;
	for(k=STARTOFFZ+1 ; k<=nnz+STARTOFFZ-2 ; k++){
		for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
			//for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
			m = nx*(ny*k + j) + i;
			if(geo_mat[m]==GEO_FLUID){
				QIN.k[QIN.numberOfBCnodes]          = kk[m];
				QIN.kN[QIN.numberOfBCnodes]         = kk[m+1];
				rhoBC[QIN.numberOfBCnodes]          = (real)0.f;
				vx[QIN.numberOfBCnodes]             = (real)0.f;
				vy[QIN.numberOfBCnodes]             = (real)0.f;
				//vz[QIN.numberOfBCnodes]             = u0;
				vz[QIN.numberOfBCnodes]             = (real)(u0*2.f)*((-4.f*i*i + nnx*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*nnx + 8.f*STARTOFFX))*(-4.f*j*j + nny*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*nny + 8.f*STARTOFFY)))/((2.f - nnx)*(2.f - nnx)*(2.f - nny)*(2.f - nny));
				//vz[QIN.numberOfBCnodes]             =  (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);
				//vz[QIN.numberOfBCnodes]             = (real)(16.f*(u0*2.f)*i*j*(nx-i)*(ny-j))/(nx*nx*ny*ny);
				deltaVz[QIN.numberOfBCnodes]        = (real)0.f;
				Q.q27[E   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[W   ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[N   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[S   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[T   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[B   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[NE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[SW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[SE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[NW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[REST][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TNE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TSW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TSE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TNW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BNE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BSW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BSE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BNW ][QIN.numberOfBCnodes] = (real)1.f;
				QIN.numberOfBCnodes++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void findKforQPressX0(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQPressX0() is deprecated! - see comment above for more information");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
	//int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
	//int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
	int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
	real ON[27];
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //unsigned int mm;
	unsigned int i, j, k, m, l;
	real test = (real) 0.f;
	int nx                        = para->getParH(lev)->nx;
	int ny                        = para->getParH(lev)->ny;
	//unsigned int nnx              = para->getParH(lev)->gridNX;
	unsigned int nny              = para->getParH(lev)->gridNY;
	unsigned int nnz              = para->getParH(lev)->gridNZ;
	int* geo_mat                  = para->getParH(lev)->geo;
	QforBoundaryConditions &QIN   = para->getParH(lev)->QpressX0;
	QIN.numberOfBCnodes = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	i=STARTOFFX+1;
	//k=nnz+STARTOFFZ-3;
	{
	for(k=STARTOFFZ+1 ; k<=nnz+STARTOFFZ-2 ; k++){
		for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
			//for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
				m = nx*(ny*k + j) + i;
				if(geo_mat[m]==GEO_FLUID){
					test = (real)0.f;
					for(l=0;l<=26;l++){
						//mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
						if(ez[l]==1){
							ON[l] = (real) 1.f;
						}
						else{
							ON[l] = (real) 0.f;
						}
						test += ON[l];
					}
					if (test>0)
					{
						QIN.numberOfBCnodes++;
					}
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
void findQPressX1(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findQPressX1() is deprecated! - see comment above for more information");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
	//int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
	//int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
	//int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int i, j, k, m;
	int nx                        = para->getParH(lev)->nx;
	int ny                        = para->getParH(lev)->ny;
	unsigned int nnx              = para->getParH(lev)->gridNX;
	unsigned int nny              = para->getParH(lev)->gridNY;
	unsigned int nnz              = para->getParH(lev)->gridNZ;
	int* geo_mat                  = para->getParH(lev)->geo;
	unsigned int* kk              = para->getParH(lev)->k;
	//unsigned int sizeQ            = para->getParH(lev)->outflowBC.numberOfBCnodes;
	unsigned int sizeQ            = para->getParH(lev)->QpressX1.numberOfBCnodes;
	real* rhoBC                = para->getParH(lev)->QpressX1.RhoBC;
	real u0                    = para->getVelocity();
	real* vx                   = para->getParH(lev)->QpressX1.Vx;
	real* vy                   = para->getParH(lev)->QpressX1.Vy;
	real* vz                   = para->getParH(lev)->QpressX1.Vz;
	real*deltaVz               = para->getParH(lev)->QpressX1.deltaVz;
	real* QQ                   = para->getParH(lev)->QpressX1.q27[0];
	QforBoundaryConditions &QIN   = para->getParH(lev)->QpressX1;
	//unsigned int nxny = nx*ny;
	QIN.numberOfBCnodes = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QforBoundaryConditions Q;
	Q.q27[E   ] = &QQ[E   *sizeQ];
	Q.q27[W   ] = &QQ[W   *sizeQ];
	Q.q27[N   ] = &QQ[N   *sizeQ];
	Q.q27[S   ] = &QQ[S   *sizeQ];
	Q.q27[T   ] = &QQ[T   *sizeQ];
	Q.q27[B   ] = &QQ[B   *sizeQ];
	Q.q27[NE  ] = &QQ[NE  *sizeQ];
	Q.q27[SW  ] = &QQ[SW  *sizeQ];
	Q.q27[SE  ] = &QQ[SE  *sizeQ];
	Q.q27[NW  ] = &QQ[NW  *sizeQ];
	Q.q27[TE  ] = &QQ[TE  *sizeQ];
	Q.q27[BW  ] = &QQ[BW  *sizeQ];
	Q.q27[BE  ] = &QQ[BE  *sizeQ];
	Q.q27[TW  ] = &QQ[TW  *sizeQ];
	Q.q27[TN  ] = &QQ[TN  *sizeQ];
	Q.q27[BS  ] = &QQ[BS  *sizeQ];
	Q.q27[BN  ] = &QQ[BN  *sizeQ];
	Q.q27[TS  ] = &QQ[TS  *sizeQ];
	Q.q27[REST] = &QQ[REST*sizeQ];
	Q.q27[TNE ] = &QQ[TNE *sizeQ];
	Q.q27[TSW ] = &QQ[TSW *sizeQ];
	Q.q27[TSE ] = &QQ[TSE *sizeQ];
	Q.q27[TNW ] = &QQ[TNW *sizeQ];
	Q.q27[BNE ] = &QQ[BNE *sizeQ];
	Q.q27[BSW ] = &QQ[BSW *sizeQ];
	Q.q27[BSE ] = &QQ[BSE *sizeQ];
	Q.q27[BNW ] = &QQ[BNW *sizeQ];


	//unsigned int li = ((nnx+STARTOFFX-2)-(STARTOFFX+1)-1);
	//unsigned int lj = ((nny+STARTOFFY-2)-(STARTOFFY+1)-1);

	i=nnx+STARTOFFX-3;
	//k=nnz+STARTOFFZ-3;
	for(k=STARTOFFZ+1 ; k<=nnz+STARTOFFZ-2 ; k++){
		for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
			//for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
			m = nx*(ny*k + j) + i;
			if(geo_mat[m]==GEO_FLUID){
				QIN.k[QIN.numberOfBCnodes]          = kk[m];
				QIN.kN[QIN.numberOfBCnodes]         = kk[m-1];
				rhoBC[QIN.numberOfBCnodes]          = (real)0.f;
				vx[QIN.numberOfBCnodes]             = (real)0.f;
				vy[QIN.numberOfBCnodes]             = (real)0.f;
				//vz[QIN.numberOfBCnodes]             = u0;
				vz[QIN.numberOfBCnodes]             = (real)(u0*2.f)*((-4.f*i*i + nnx*(-2.f - 4.f*STARTOFFX) - 4.f*(-1.5f + STARTOFFX)*(0.5f + STARTOFFX) + i*(-4.f + 4.f*nnx + 8.f*STARTOFFX))*(-4.f*j*j + nny*(-2.f - 4.f*STARTOFFY) - 4.f*(-1.5f + STARTOFFY)*(0.5f + STARTOFFY) + j*(-4.f + 4.f*nny + 8.f*STARTOFFY)))/((2.f - nnx)*(2.f - nnx)*(2.f - nny)*(2.f - nny));
				//vz[QIN.numberOfBCnodes]             =  (real)(16.f*(u0*2.f)*(i-(STARTOFFX+1)-0.5f)*(li-1.5f-(i-(STARTOFFX+1)))*(j-(STARTOFFY+1)-0.5f)*(lj-1.5f-(j-(STARTOFFY+1))))/(li*lj*li*lj);
				//vz[QIN.numberOfBCnodes]             = (real)(16.f*(u0*2.f)*i*j*(nx-i)*(ny-j))/(nx*nx*ny*ny);
				deltaVz[QIN.numberOfBCnodes]        = (real)0.f;
				Q.q27[E   ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[W   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[N   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[S   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[T   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[B   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[NE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[SW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[SE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[NW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[REST][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TNE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TSW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[TSE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[TNW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BNE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BSW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[BSE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[BNW ][QIN.numberOfBCnodes] = (real)-1.f;
				QIN.numberOfBCnodes++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void findKforQPressX1(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQPressX1() is deprecated! - see comment above for more information");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////  E   W   N   S   T   B  NE  SW  SE  NW  TE  BW  BE  TW  TN  BS  BN  TS ZERO TNE BNE TSE BSE TNW BNW TSW BSW  ////////////////////////
	//int   ex[27]={  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,   0,  1,  1,  1,  1, -1, -1, -1, -1};
	//int   ey[27]={  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,   0,  1,  1, -1, -1,  1,  1, -1, -1};
	int   ez[27]={  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,   0,  1, -1,  1, -1,  1, -1,  1, -1};
	real ON[27];
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //unsigned int mm;
	unsigned int i, j, k, m, l;
	real test = (real) 0.f;
	int nx                        = para->getParH(lev)->nx;
	int ny                        = para->getParH(lev)->ny;
	unsigned int nnx              = para->getParH(lev)->gridNX;
	unsigned int nny              = para->getParH(lev)->gridNY;
	unsigned int nnz              = para->getParH(lev)->gridNZ;
	int* geo_mat                  = para->getParH(lev)->geo;
	QforBoundaryConditions &QIN   = para->getParH(lev)->QpressX1;
	QIN.numberOfBCnodes = 0;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	i=nnx+STARTOFFX-3;
	//k=nnz+STARTOFFZ-3;
	{
		for(k=STARTOFFZ+1 ; k<=nnz+STARTOFFZ-2 ; k++){
			for(j=STARTOFFY+1 ; j<=nny+STARTOFFY-2 ; j++){
				//for(i=STARTOFFX+1; i<=nnx+STARTOFFX-2 ; i++){
				m = nx*(ny*k + j) + i;
				if(geo_mat[m]==GEO_FLUID){
					test = (real)0.f;
					for(l=0;l<=26;l++){
						//mm = nx*(ny*(k+ez[l]) + (j+ey[l])) + (i+ex[l]);
						if(ez[l]==1){
							ON[l] = (real) 1.f;
						}
						else{
							ON[l] = (real) 0.f;
						}
						test += ON[l];
					}
					if (test>0)
					{
						QIN.numberOfBCnodes++;
					}
				}
			}
		}
	}
}
