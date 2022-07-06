#include "FindQ/FindQ.h"
#include <stdexcept>


////////////////////////////////////////////////////////////////////////////////
void findQ(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

   throw std::runtime_error("findQ() is deprecated! - see comment above for more information");

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
   Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
   Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
   Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
   Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
   Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
   Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
   Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
   Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
   Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
   Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
   Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
   Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
   Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
   Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
   Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
   Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
   Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
   Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
   Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
   Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
   Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
   Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
   Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
   Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
   Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
   Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
   Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
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

                  //Q.q27[dirE   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirW   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirN   ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirS   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirT   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirB   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirNE  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirSW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirSE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirNW  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirTE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirBW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirBE  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirTW  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirTN  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirBS  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirBN  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirTS  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirZERO][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirTNE ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirTSW ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirTSE ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirTNW ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirBNE ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dirBSW ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirBSE ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dirBNW ][QIN.numberOfBCnodes] = 0.f;

                  //Q.q27[dirE   ][QIN.numberOfBCnodes] = ON[dirW   ];
                  //Q.q27[dirW   ][QIN.numberOfBCnodes] = ON[dirE   ];
                  //Q.q27[dirN   ][QIN.numberOfBCnodes] = ON[dirS   ];
                  //Q.q27[dirS   ][QIN.numberOfBCnodes] = ON[dirN   ];
                  //Q.q27[dirT   ][QIN.numberOfBCnodes] = ON[dirB   ];
                  //Q.q27[dirB   ][QIN.numberOfBCnodes] = ON[dirT   ];
                  //Q.q27[dirNE  ][QIN.numberOfBCnodes] = ON[dirSW  ];
                  //Q.q27[dirSW  ][QIN.numberOfBCnodes] = ON[dirNE  ];
                  //Q.q27[dirSE  ][QIN.numberOfBCnodes] = ON[dirNW  ];
                  //Q.q27[dirNW  ][QIN.numberOfBCnodes] = ON[dirSE  ];
                  //Q.q27[dirTE  ][QIN.numberOfBCnodes] = ON[dirBW  ];
                  //Q.q27[dirBW  ][QIN.numberOfBCnodes] = ON[dirTE  ];
                  //Q.q27[dirBE  ][QIN.numberOfBCnodes] = ON[dirTW  ];
                  //Q.q27[dirTW  ][QIN.numberOfBCnodes] = ON[dirBE  ];
                  //Q.q27[dirTN  ][QIN.numberOfBCnodes] = ON[dirBS  ];
                  //Q.q27[dirBS  ][QIN.numberOfBCnodes] = ON[dirTN  ];
                  //Q.q27[dirBN  ][QIN.numberOfBCnodes] = ON[dirTS  ];
                  //Q.q27[dirTS  ][QIN.numberOfBCnodes] = ON[dirBN  ];
                  //Q.q27[dirZERO][QIN.numberOfBCnodes] = ON[dirZERO];
                  //Q.q27[dirTNE ][QIN.numberOfBCnodes] = ON[dirBSW ];
                  //Q.q27[dirTSW ][QIN.numberOfBCnodes] = ON[dirBNE ];
                  //Q.q27[dirTSE ][QIN.numberOfBCnodes] = ON[dirBNW ];
                  //Q.q27[dirTNW ][QIN.numberOfBCnodes] = ON[dirBSE ];
                  //Q.q27[dirBNE ][QIN.numberOfBCnodes] = ON[dirTSW ];
                  //Q.q27[dirBSW ][QIN.numberOfBCnodes] = ON[dirTNE ];
                  //Q.q27[dirBSE ][QIN.numberOfBCnodes] = ON[dirTNW ];
                  //Q.q27[dirBNW ][QIN.numberOfBCnodes] = ON[dirTSE ];

                  Q.q27[dirE   ][QIN.numberOfBCnodes] = ON[dirE   ];
                  Q.q27[dirW   ][QIN.numberOfBCnodes] = ON[dirW   ];
                  Q.q27[dirN   ][QIN.numberOfBCnodes] = ON[dirN   ];
                  Q.q27[dirS   ][QIN.numberOfBCnodes] = ON[dirS   ];
                  Q.q27[dirT   ][QIN.numberOfBCnodes] = ON[dirT   ];
                  Q.q27[dirB   ][QIN.numberOfBCnodes] = ON[dirB   ];
                  Q.q27[dirNE  ][QIN.numberOfBCnodes] = ON[dirNE  ];
                  Q.q27[dirSW  ][QIN.numberOfBCnodes] = ON[dirSW  ];
                  Q.q27[dirSE  ][QIN.numberOfBCnodes] = ON[dirSE  ];
                  Q.q27[dirNW  ][QIN.numberOfBCnodes] = ON[dirNW  ];
                  Q.q27[dirTE  ][QIN.numberOfBCnodes] = ON[dirTE  ];
                  Q.q27[dirBW  ][QIN.numberOfBCnodes] = ON[dirBW  ];
                  Q.q27[dirBE  ][QIN.numberOfBCnodes] = ON[dirBE  ];
                  Q.q27[dirTW  ][QIN.numberOfBCnodes] = ON[dirTW  ];
                  Q.q27[dirTN  ][QIN.numberOfBCnodes] = ON[dirTN  ];
                  Q.q27[dirBS  ][QIN.numberOfBCnodes] = ON[dirBS  ];
                  Q.q27[dirBN  ][QIN.numberOfBCnodes] = ON[dirBN  ];
                  Q.q27[dirTS  ][QIN.numberOfBCnodes] = ON[dirTS  ];
                  Q.q27[dirZERO][QIN.numberOfBCnodes] = ON[dirZERO];
                  Q.q27[dirTNE ][QIN.numberOfBCnodes] = ON[dirTNE ];
                  Q.q27[dirTSW ][QIN.numberOfBCnodes] = ON[dirTSW ];
                  Q.q27[dirTSE ][QIN.numberOfBCnodes] = ON[dirTSE ];
                  Q.q27[dirTNW ][QIN.numberOfBCnodes] = ON[dirTNW ];
                  Q.q27[dirBNE ][QIN.numberOfBCnodes] = ON[dirBNE ];
                  Q.q27[dirBSW ][QIN.numberOfBCnodes] = ON[dirBSW ];
                  Q.q27[dirBSE ][QIN.numberOfBCnodes] = ON[dirBSE ];
                  Q.q27[dirBNW ][QIN.numberOfBCnodes] = ON[dirBNW ];

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

   throw std::runtime_error("findKforQ() is deprecated! - see comment above for more information");

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
   Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
   Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
   Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
   Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
   Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
   Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
   Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
   Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
   Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
   Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
   Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
   Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
   Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
   Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
   Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
   Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
   Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
   Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
   Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
   Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
   Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
   Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
   Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
   Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
   Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
   Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
   Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];

   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

   throw std::runtime_error("findQ_MG() is deprecated! - see comment above for more information");

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

                  Q.q27[dirE   ][QIN.numberOfBCnodes] = ON[dirE   ];
                  Q.q27[dirW   ][QIN.numberOfBCnodes] = ON[dirW   ];
                  Q.q27[dirN   ][QIN.numberOfBCnodes] = ON[dirN   ];
                  Q.q27[dirS   ][QIN.numberOfBCnodes] = ON[dirS   ];
                  Q.q27[dirT   ][QIN.numberOfBCnodes] = ON[dirT   ];
                  Q.q27[dirB   ][QIN.numberOfBCnodes] = ON[dirB   ];
                  Q.q27[dirNE  ][QIN.numberOfBCnodes] = ON[dirNE  ];
                  Q.q27[dirSW  ][QIN.numberOfBCnodes] = ON[dirSW  ];
                  Q.q27[dirSE  ][QIN.numberOfBCnodes] = ON[dirSE  ];
                  Q.q27[dirNW  ][QIN.numberOfBCnodes] = ON[dirNW  ];
                  Q.q27[dirTE  ][QIN.numberOfBCnodes] = ON[dirTE  ];
                  Q.q27[dirBW  ][QIN.numberOfBCnodes] = ON[dirBW  ];
                  Q.q27[dirBE  ][QIN.numberOfBCnodes] = ON[dirBE  ];
                  Q.q27[dirTW  ][QIN.numberOfBCnodes] = ON[dirTW  ];
                  Q.q27[dirTN  ][QIN.numberOfBCnodes] = ON[dirTN  ];
                  Q.q27[dirBS  ][QIN.numberOfBCnodes] = ON[dirBS  ];
                  Q.q27[dirBN  ][QIN.numberOfBCnodes] = ON[dirBN  ];
                  Q.q27[dirTS  ][QIN.numberOfBCnodes] = ON[dirTS  ];
                  Q.q27[dirZERO][QIN.numberOfBCnodes] = ON[dirZERO];
                  Q.q27[dirTNE ][QIN.numberOfBCnodes] = ON[dirTNE ];
                  Q.q27[dirTSW ][QIN.numberOfBCnodes] = ON[dirTSW ];
                  Q.q27[dirTSE ][QIN.numberOfBCnodes] = ON[dirTSE ];
                  Q.q27[dirTNW ][QIN.numberOfBCnodes] = ON[dirTNW ];
                  Q.q27[dirBNE ][QIN.numberOfBCnodes] = ON[dirBNE ];
                  Q.q27[dirBSW ][QIN.numberOfBCnodes] = ON[dirBSW ];
                  Q.q27[dirBSE ][QIN.numberOfBCnodes] = ON[dirBSE ];
                  Q.q27[dirBNW ][QIN.numberOfBCnodes] = ON[dirBNW ];

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

   throw std::runtime_error("findKforQ_MG() is deprecated! - see comment above for more information");

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
   throw std::runtime_error("findQInflow() is deprecated! - see comment above for more information");

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
   Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
   Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
   Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
   Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
   Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
   Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
   Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
   Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
   Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
   Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
   Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
   Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
   Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
   Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
   Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
   Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
   Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
   Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
   Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
   Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
   Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
   Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
   Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
   Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
   Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
   Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
   Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];
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
               //Q.q27[dirE   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirW   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirN   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirS   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirT   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirB   ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirNE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirSW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirSE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirNW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTE  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirBW  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirBE  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirTW  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTN  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirBS  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirBN  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirTS  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirZERO][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTNE ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTSW ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTSE ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirTNW ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dirBNE ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirBSW ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirBSE ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dirBNW ][QIN.numberOfBCnodes] = (real)1.f;
			   //////////////////////////////////////////////////////////////////////////


               // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

			   Q.q27[dirE   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirW   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirN   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirS   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirT   ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirB   ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirNE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirSW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirSE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirNW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirTE  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirBW  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirBE  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirTW  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirTN  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirBS  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirBN  ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirTS  ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirZERO][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirTNE ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirTSW ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirTSE ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirTNW ][QIN.numberOfBCnodes] = (real)1.f;
			   Q.q27[dirBNE ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirBSW ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirBSE ][QIN.numberOfBCnodes] = (real)-1.f;
			   Q.q27[dirBNW ][QIN.numberOfBCnodes] = (real)-1.f;
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

   //               Q.q27[dirE   ][QIN.numberOfBCnodes] = ON[dirE   ];
   //               Q.q27[dirW   ][QIN.numberOfBCnodes] = ON[dirW   ];
   //               Q.q27[dirN   ][QIN.numberOfBCnodes] = ON[dirN   ];
   //               Q.q27[dirS   ][QIN.numberOfBCnodes] = ON[dirS   ];
   //               Q.q27[dirT   ][QIN.numberOfBCnodes] = ON[dirT   ];
   //               Q.q27[dirB   ][QIN.numberOfBCnodes] = ON[dirB   ];
   //               Q.q27[dirNE  ][QIN.numberOfBCnodes] = ON[dirNE  ];
   //               Q.q27[dirSW  ][QIN.numberOfBCnodes] = ON[dirSW  ];
   //               Q.q27[dirSE  ][QIN.numberOfBCnodes] = ON[dirSE  ];
   //               Q.q27[dirNW  ][QIN.numberOfBCnodes] = ON[dirNW  ];
   //               Q.q27[dirTE  ][QIN.numberOfBCnodes] = ON[dirTE  ];
   //               Q.q27[dirBW  ][QIN.numberOfBCnodes] = ON[dirBW  ];
   //               Q.q27[dirBE  ][QIN.numberOfBCnodes] = ON[dirBE  ];
   //               Q.q27[dirTW  ][QIN.numberOfBCnodes] = ON[dirTW  ];
   //               Q.q27[dirTN  ][QIN.numberOfBCnodes] = ON[dirTN  ];
   //               Q.q27[dirBS  ][QIN.numberOfBCnodes] = ON[dirBS  ];
   //               Q.q27[dirBN  ][QIN.numberOfBCnodes] = ON[dirBN  ];
   //               Q.q27[dirTS  ][QIN.numberOfBCnodes] = ON[dirTS  ];
   //               Q.q27[dirZERO][QIN.numberOfBCnodes] = ON[dirZERO];
   //               Q.q27[dirTNE ][QIN.numberOfBCnodes] = ON[dirTNE ];
   //               Q.q27[dirTSW ][QIN.numberOfBCnodes] = ON[dirTSW ];
   //               Q.q27[dirTSE ][QIN.numberOfBCnodes] = ON[dirTSE ];
   //               Q.q27[dirTNW ][QIN.numberOfBCnodes] = ON[dirTNW ];
   //               Q.q27[dirBNE ][QIN.numberOfBCnodes] = ON[dirBNE ];
   //               Q.q27[dirBSW ][QIN.numberOfBCnodes] = ON[dirBSW ];
   //               Q.q27[dirBSE ][QIN.numberOfBCnodes] = ON[dirBSE ];
   //               Q.q27[dirBNW ][QIN.numberOfBCnodes] = ON[dirBNW ];

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

   //               Q.q27[dirE   ][QIN.numberOfBCnodes] = ON[dirE   ];
   //               Q.q27[dirW   ][QIN.numberOfBCnodes] = ON[dirW   ];
   //               Q.q27[dirN   ][QIN.numberOfBCnodes] = ON[dirN   ];
   //               Q.q27[dirS   ][QIN.numberOfBCnodes] = ON[dirS   ];
   //               Q.q27[dirT   ][QIN.numberOfBCnodes] = ON[dirT   ];
   //               Q.q27[dirB   ][QIN.numberOfBCnodes] = ON[dirB   ];
   //               Q.q27[dirNE  ][QIN.numberOfBCnodes] = ON[dirNE  ];
   //               Q.q27[dirSW  ][QIN.numberOfBCnodes] = ON[dirSW  ];
   //               Q.q27[dirSE  ][QIN.numberOfBCnodes] = ON[dirSE  ];
   //               Q.q27[dirNW  ][QIN.numberOfBCnodes] = ON[dirNW  ];
   //               Q.q27[dirTE  ][QIN.numberOfBCnodes] = ON[dirTE  ];
   //               Q.q27[dirBW  ][QIN.numberOfBCnodes] = ON[dirBW  ];
   //               Q.q27[dirBE  ][QIN.numberOfBCnodes] = ON[dirBE  ];
   //               Q.q27[dirTW  ][QIN.numberOfBCnodes] = ON[dirTW  ];
   //               Q.q27[dirTN  ][QIN.numberOfBCnodes] = ON[dirTN  ];
   //               Q.q27[dirBS  ][QIN.numberOfBCnodes] = ON[dirBS  ];
   //               Q.q27[dirBN  ][QIN.numberOfBCnodes] = ON[dirBN  ];
   //               Q.q27[dirTS  ][QIN.numberOfBCnodes] = ON[dirTS  ];
   //               Q.q27[dirZERO][QIN.numberOfBCnodes] = ON[dirZERO];
   //               Q.q27[dirTNE ][QIN.numberOfBCnodes] = ON[dirTNE ];
   //               Q.q27[dirTSW ][QIN.numberOfBCnodes] = ON[dirTSW ];
   //               Q.q27[dirTSE ][QIN.numberOfBCnodes] = ON[dirTSE ];
   //               Q.q27[dirTNW ][QIN.numberOfBCnodes] = ON[dirTNW ];
   //               Q.q27[dirBNE ][QIN.numberOfBCnodes] = ON[dirBNE ];
   //               Q.q27[dirBSW ][QIN.numberOfBCnodes] = ON[dirBSW ];
   //               Q.q27[dirBSE ][QIN.numberOfBCnodes] = ON[dirBSE ];
   //               Q.q27[dirBNW ][QIN.numberOfBCnodes] = ON[dirBNW ];

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
   throw std::runtime_error("findKforQInflow() is deprecated! - see comment above for more information");

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
   throw std::runtime_error("findQOutflow() is deprecated! - see comment above for more information");

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
   Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
   Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
   Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
   Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
   Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
   Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
   Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
   Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
   Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
   Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
   Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
   Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
   Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
   Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
   Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
   Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
   Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
   Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
   Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
   Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
   Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
   Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
   Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
   Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
   Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
   Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
   Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];


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
               Q.q27[dirE   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirW   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirN   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirS   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirT   ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirB   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirNE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirSW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirSE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirNW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirTE  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirBW  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirBE  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirTW  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirTN  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirBS  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirBN  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirTS  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirZERO][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirTNE ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirTSW ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirTSE ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirTNW ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dirBNE ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirBSW ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirBSE ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dirBNW ][QIN.numberOfBCnodes] = (real)-1.f;
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
   throw std::runtime_error("findKforQOutflow() is deprecated! - see comment above for more information");

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
//   QN.q27[dirE   ] = &QQN[dirE   *sizeQN];
//   QN.q27[dirW   ] = &QQN[dirW   *sizeQN];
//   QN.q27[dirN   ] = &QQN[dirN   *sizeQN];
//   QN.q27[dirS   ] = &QQN[dirS   *sizeQN];
//   QN.q27[dirT   ] = &QQN[dirT   *sizeQN];
//   QN.q27[dirB   ] = &QQN[dirB   *sizeQN];
//   QN.q27[dirNE  ] = &QQN[dirNE  *sizeQN];
//   QN.q27[dirSW  ] = &QQN[dirSW  *sizeQN];
//   QN.q27[dirSE  ] = &QQN[dirSE  *sizeQN];
//   QN.q27[dirNW  ] = &QQN[dirNW  *sizeQN];
//   QN.q27[dirTE  ] = &QQN[dirTE  *sizeQN];
//   QN.q27[dirBW  ] = &QQN[dirBW  *sizeQN];
//   QN.q27[dirBE  ] = &QQN[dirBE  *sizeQN];
//   QN.q27[dirTW  ] = &QQN[dirTW  *sizeQN];
//   QN.q27[dirTN  ] = &QQN[dirTN  *sizeQN];
//   QN.q27[dirBS  ] = &QQN[dirBS  *sizeQN];
//   QN.q27[dirBN  ] = &QQN[dirBN  *sizeQN];
//   QN.q27[dirTS  ] = &QQN[dirTS  *sizeQN];
//   QN.q27[dirZERO] = &QQN[dirZERO*sizeQN];
//   QN.q27[dirTNE ] = &QQN[dirTNE *sizeQN];
//   QN.q27[dirTSW ] = &QQN[dirTSW *sizeQN];
//   QN.q27[dirTSE ] = &QQN[dirTSE *sizeQN];
//   QN.q27[dirTNW ] = &QQN[dirTNW *sizeQN];
//   QN.q27[dirBNE ] = &QQN[dirBNE *sizeQN];
//   QN.q27[dirBSW ] = &QQN[dirBSW *sizeQN];
//   QN.q27[dirBSE ] = &QQN[dirBSE *sizeQN];
//   QN.q27[dirBNW ] = &QQN[dirBNW *sizeQN];
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
   throw std::runtime_error("findKforQPressX0() is deprecated! - see comment above for more information");

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
	Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
	Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
	Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
	Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
	Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
	Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
	Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
	Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
	Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
	Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
	Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
	Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
	Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
	Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
	Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
	Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
	Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
	Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
	Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
	Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
	Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
	Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
	Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
	Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
	Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
	Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
	Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];


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
				Q.q27[dirE   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirW   ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirN   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirS   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirT   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirB   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirNE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirSW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirSE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirNW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBE  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTW  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirZERO][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTNE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTSW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTSE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTNW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBNE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBSW ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBSE ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBNW ][QIN.numberOfBCnodes] = (real)1.f;
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
   throw std::runtime_error("findKforQPressX0() is deprecated! - see comment above for more information");

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
   throw std::runtime_error("findQPressX1() is deprecated! - see comment above for more information");

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
	Q.q27[dirE   ] = &QQ[dirE   *sizeQ];
	Q.q27[dirW   ] = &QQ[dirW   *sizeQ];
	Q.q27[dirN   ] = &QQ[dirN   *sizeQ];
	Q.q27[dirS   ] = &QQ[dirS   *sizeQ];
	Q.q27[dirT   ] = &QQ[dirT   *sizeQ];
	Q.q27[dirB   ] = &QQ[dirB   *sizeQ];
	Q.q27[dirNE  ] = &QQ[dirNE  *sizeQ];
	Q.q27[dirSW  ] = &QQ[dirSW  *sizeQ];
	Q.q27[dirSE  ] = &QQ[dirSE  *sizeQ];
	Q.q27[dirNW  ] = &QQ[dirNW  *sizeQ];
	Q.q27[dirTE  ] = &QQ[dirTE  *sizeQ];
	Q.q27[dirBW  ] = &QQ[dirBW  *sizeQ];
	Q.q27[dirBE  ] = &QQ[dirBE  *sizeQ];
	Q.q27[dirTW  ] = &QQ[dirTW  *sizeQ];
	Q.q27[dirTN  ] = &QQ[dirTN  *sizeQ];
	Q.q27[dirBS  ] = &QQ[dirBS  *sizeQ];
	Q.q27[dirBN  ] = &QQ[dirBN  *sizeQ];
	Q.q27[dirTS  ] = &QQ[dirTS  *sizeQ];
	Q.q27[dirZERO] = &QQ[dirZERO*sizeQ];
	Q.q27[dirTNE ] = &QQ[dirTNE *sizeQ];
	Q.q27[dirTSW ] = &QQ[dirTSW *sizeQ];
	Q.q27[dirTSE ] = &QQ[dirTSE *sizeQ];
	Q.q27[dirTNW ] = &QQ[dirTNW *sizeQ];
	Q.q27[dirBNE ] = &QQ[dirBNE *sizeQ];
	Q.q27[dirBSW ] = &QQ[dirBSW *sizeQ];
	Q.q27[dirBSE ] = &QQ[dirBSE *sizeQ];
	Q.q27[dirBNW ] = &QQ[dirBNW *sizeQ];


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
				Q.q27[dirE   ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirW   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirN   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirS   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirT   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirB   ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirNE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirSW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirSE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirNW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBE  ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTW  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBN  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTS  ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirZERO][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTNE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTSW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirTSE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirTNW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBNE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBSW ][QIN.numberOfBCnodes] = (real)-1.f;
				Q.q27[dirBSE ][QIN.numberOfBCnodes] = (real)1.f;
				Q.q27[dirBNW ][QIN.numberOfBCnodes] = (real)-1.f;
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
   throw std::runtime_error("findKforQPressX1() is deprecated! - see comment above for more information");

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
