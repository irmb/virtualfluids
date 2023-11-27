#include "FindQ/FindQ.h"
#include <logger/Logger.h>
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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
   Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
   Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
   Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
   Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
   Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
   Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
   Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
   Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
   Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
   Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
   Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
   Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
   Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
   Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
   Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
   Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
   Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
   Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
   Q.q27[d000] = &QQ[d000*sizeQ];
   Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
   Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
   Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
   Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
   Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
   Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
   Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
   Q.q27[dMPM ] = &QQ[dMPM *sizeQ];
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

                  //Q.q27[dP00   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dM00   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d0P0   ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[d0M0   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d00P   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d00M   ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dPP0  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dMM0  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dPM0  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dMP0  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dP0P  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dM0M  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dP0M  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dM0P  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d0PP  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[d0MM  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d0PM  ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[d0MP  ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[d000][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dPPP ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dMMP ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dPMP ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dMPP ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dPPM ][QIN.numberOfBCnodes] = 0.f;
                  //Q.q27[dMMM ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dPMM ][QIN.numberOfBCnodes] = -1.f;
                  //Q.q27[dMPM ][QIN.numberOfBCnodes] = 0.f;

                  //Q.q27[dP00   ][QIN.numberOfBCnodes] = ON[dM00   ];
                  //Q.q27[dM00   ][QIN.numberOfBCnodes] = ON[dP00   ];
                  //Q.q27[d0P0   ][QIN.numberOfBCnodes] = ON[d0M0   ];
                  //Q.q27[d0M0   ][QIN.numberOfBCnodes] = ON[d0P0   ];
                  //Q.q27[d00P   ][QIN.numberOfBCnodes] = ON[d00M   ];
                  //Q.q27[d00M   ][QIN.numberOfBCnodes] = ON[d00P   ];
                  //Q.q27[dPP0  ][QIN.numberOfBCnodes] = ON[dMM0  ];
                  //Q.q27[dMM0  ][QIN.numberOfBCnodes] = ON[dPP0  ];
                  //Q.q27[dPM0  ][QIN.numberOfBCnodes] = ON[dMP0  ];
                  //Q.q27[dMP0  ][QIN.numberOfBCnodes] = ON[dPM0  ];
                  //Q.q27[dP0P  ][QIN.numberOfBCnodes] = ON[dM0M  ];
                  //Q.q27[dM0M  ][QIN.numberOfBCnodes] = ON[dP0P  ];
                  //Q.q27[dP0M  ][QIN.numberOfBCnodes] = ON[dM0P  ];
                  //Q.q27[dM0P  ][QIN.numberOfBCnodes] = ON[dP0M  ];
                  //Q.q27[d0PP  ][QIN.numberOfBCnodes] = ON[d0MM  ];
                  //Q.q27[d0MM  ][QIN.numberOfBCnodes] = ON[d0PP  ];
                  //Q.q27[d0PM  ][QIN.numberOfBCnodes] = ON[d0MP  ];
                  //Q.q27[d0MP  ][QIN.numberOfBCnodes] = ON[d0PM  ];
                  //Q.q27[d000][QIN.numberOfBCnodes] = ON[d000];
                  //Q.q27[dPPP ][QIN.numberOfBCnodes] = ON[dMMM ];
                  //Q.q27[dMMP ][QIN.numberOfBCnodes] = ON[dPPM ];
                  //Q.q27[dPMP ][QIN.numberOfBCnodes] = ON[dMPM ];
                  //Q.q27[dMPP ][QIN.numberOfBCnodes] = ON[dPMM ];
                  //Q.q27[dPPM ][QIN.numberOfBCnodes] = ON[dMMP ];
                  //Q.q27[dMMM ][QIN.numberOfBCnodes] = ON[dPPP ];
                  //Q.q27[dPMM ][QIN.numberOfBCnodes] = ON[dMPP ];
                  //Q.q27[dMPM ][QIN.numberOfBCnodes] = ON[dPMP ];

                  Q.q27[dP00   ][QIN.numberOfBCnodes] = ON[dP00   ];
                  Q.q27[dM00   ][QIN.numberOfBCnodes] = ON[dM00   ];
                  Q.q27[d0P0   ][QIN.numberOfBCnodes] = ON[d0P0   ];
                  Q.q27[d0M0   ][QIN.numberOfBCnodes] = ON[d0M0   ];
                  Q.q27[d00P   ][QIN.numberOfBCnodes] = ON[d00P   ];
                  Q.q27[d00M   ][QIN.numberOfBCnodes] = ON[d00M   ];
                  Q.q27[dPP0  ][QIN.numberOfBCnodes] = ON[dPP0  ];
                  Q.q27[dMM0  ][QIN.numberOfBCnodes] = ON[dMM0  ];
                  Q.q27[dPM0  ][QIN.numberOfBCnodes] = ON[dPM0  ];
                  Q.q27[dMP0  ][QIN.numberOfBCnodes] = ON[dMP0  ];
                  Q.q27[dP0P  ][QIN.numberOfBCnodes] = ON[dP0P  ];
                  Q.q27[dM0M  ][QIN.numberOfBCnodes] = ON[dM0M  ];
                  Q.q27[dP0M  ][QIN.numberOfBCnodes] = ON[dP0M  ];
                  Q.q27[dM0P  ][QIN.numberOfBCnodes] = ON[dM0P  ];
                  Q.q27[d0PP  ][QIN.numberOfBCnodes] = ON[d0PP  ];
                  Q.q27[d0MM  ][QIN.numberOfBCnodes] = ON[d0MM  ];
                  Q.q27[d0PM  ][QIN.numberOfBCnodes] = ON[d0PM  ];
                  Q.q27[d0MP  ][QIN.numberOfBCnodes] = ON[d0MP  ];
                  Q.q27[d000][QIN.numberOfBCnodes] = ON[d000];
                  Q.q27[dPPP ][QIN.numberOfBCnodes] = ON[dPPP ];
                  Q.q27[dMMP ][QIN.numberOfBCnodes] = ON[dMMP ];
                  Q.q27[dPMP ][QIN.numberOfBCnodes] = ON[dPMP ];
                  Q.q27[dMPP ][QIN.numberOfBCnodes] = ON[dMPP ];
                  Q.q27[dPPM ][QIN.numberOfBCnodes] = ON[dPPM ];
                  Q.q27[dMMM ][QIN.numberOfBCnodes] = ON[dMMM ];
                  Q.q27[dPMM ][QIN.numberOfBCnodes] = ON[dPMM ];
                  Q.q27[dMPM ][QIN.numberOfBCnodes] = ON[dMPM ];

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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
   Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
   Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
   Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
   Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
   Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
   Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
   Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
   Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
   Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
   Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
   Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
   Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
   Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
   Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
   Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
   Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
   Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
   Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
   Q.q27[d000] = &QQ[d000*sizeQ];
   Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
   Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
   Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
   Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
   Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
   Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
   Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
   Q.q27[dMPM ] = &QQ[dMPM *sizeQ];

   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

    VF_LOG_CRITICAL("findQ_MG() is deprecated! - see comment above for more information");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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

                  Q.q27[dP00   ][QIN.numberOfBCnodes] = ON[dP00   ];
                  Q.q27[dM00   ][QIN.numberOfBCnodes] = ON[dM00   ];
                  Q.q27[d0P0   ][QIN.numberOfBCnodes] = ON[d0P0   ];
                  Q.q27[d0M0   ][QIN.numberOfBCnodes] = ON[d0M0   ];
                  Q.q27[d00P   ][QIN.numberOfBCnodes] = ON[d00P   ];
                  Q.q27[d00M   ][QIN.numberOfBCnodes] = ON[d00M   ];
                  Q.q27[dPP0  ][QIN.numberOfBCnodes] = ON[dPP0  ];
                  Q.q27[dMM0  ][QIN.numberOfBCnodes] = ON[dMM0  ];
                  Q.q27[dPM0  ][QIN.numberOfBCnodes] = ON[dPM0  ];
                  Q.q27[dMP0  ][QIN.numberOfBCnodes] = ON[dMP0  ];
                  Q.q27[dP0P  ][QIN.numberOfBCnodes] = ON[dP0P  ];
                  Q.q27[dM0M  ][QIN.numberOfBCnodes] = ON[dM0M  ];
                  Q.q27[dP0M  ][QIN.numberOfBCnodes] = ON[dP0M  ];
                  Q.q27[dM0P  ][QIN.numberOfBCnodes] = ON[dM0P  ];
                  Q.q27[d0PP  ][QIN.numberOfBCnodes] = ON[d0PP  ];
                  Q.q27[d0MM  ][QIN.numberOfBCnodes] = ON[d0MM  ];
                  Q.q27[d0PM  ][QIN.numberOfBCnodes] = ON[d0PM  ];
                  Q.q27[d0MP  ][QIN.numberOfBCnodes] = ON[d0MP  ];
                  Q.q27[d000][QIN.numberOfBCnodes] = ON[d000];
                  Q.q27[dPPP ][QIN.numberOfBCnodes] = ON[dPPP ];
                  Q.q27[dMMP ][QIN.numberOfBCnodes] = ON[dMMP ];
                  Q.q27[dPMP ][QIN.numberOfBCnodes] = ON[dPMP ];
                  Q.q27[dMPP ][QIN.numberOfBCnodes] = ON[dMPP ];
                  Q.q27[dPPM ][QIN.numberOfBCnodes] = ON[dPPM ];
                  Q.q27[dMMM ][QIN.numberOfBCnodes] = ON[dMMM ];
                  Q.q27[dPMM ][QIN.numberOfBCnodes] = ON[dPMM ];
                  Q.q27[dMPM ][QIN.numberOfBCnodes] = ON[dMPM ];

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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
   Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
   Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
   Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
   Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
   Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
   Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
   Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
   Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
   Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
   Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
   Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
   Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
   Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
   Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
   Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
   Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
   Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
   Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
   Q.q27[d000] = &QQ[d000*sizeQ];
   Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
   Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
   Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
   Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
   Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
   Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
   Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
   Q.q27[dMPM ] = &QQ[dMPM *sizeQ];
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
               //Q.q27[dP00   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dM00   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d0P0   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d0M0   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d00P   ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d00M   ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dPP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dMM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dPM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dMP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dP0P  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dM0M  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dP0M  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dM0P  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d0PP  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d0MM  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[d0PM  ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[d0MP  ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[d000][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dPPP ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dMMP ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dPMP ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dMPP ][QIN.numberOfBCnodes] = (real)-1.f;
               //Q.q27[dPPM ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dMMM ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dPMM ][QIN.numberOfBCnodes] = (real)1.f;
               //Q.q27[dMPM ][QIN.numberOfBCnodes] = (real)1.f;
               //////////////////////////////////////////////////////////////////////////


               // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14

               Q.q27[dP00   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dM00   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0P0   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0M0   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d00P   ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d00M   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dP0P  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dM0M  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dP0M  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dM0P  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d0PP  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d0MM  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0PM  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0MP  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d000][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPPP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dMMP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dPMP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dMPP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dPPM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMMM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPMM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMPM ][QIN.numberOfBCnodes] = (real)-1.f;
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

   //               Q.q27[dP00   ][QIN.numberOfBCnodes] = ON[dP00   ];
   //               Q.q27[dM00   ][QIN.numberOfBCnodes] = ON[dM00   ];
   //               Q.q27[d0P0   ][QIN.numberOfBCnodes] = ON[d0P0   ];
   //               Q.q27[d0M0   ][QIN.numberOfBCnodes] = ON[d0M0   ];
   //               Q.q27[d00P   ][QIN.numberOfBCnodes] = ON[d00P   ];
   //               Q.q27[d00M   ][QIN.numberOfBCnodes] = ON[d00M   ];
   //               Q.q27[dPP0  ][QIN.numberOfBCnodes] = ON[dPP0  ];
   //               Q.q27[dMM0  ][QIN.numberOfBCnodes] = ON[dMM0  ];
   //               Q.q27[dPM0  ][QIN.numberOfBCnodes] = ON[dPM0  ];
   //               Q.q27[dMP0  ][QIN.numberOfBCnodes] = ON[dMP0  ];
   //               Q.q27[dP0P  ][QIN.numberOfBCnodes] = ON[dP0P  ];
   //               Q.q27[dM0M  ][QIN.numberOfBCnodes] = ON[dM0M  ];
   //               Q.q27[dP0M  ][QIN.numberOfBCnodes] = ON[dP0M  ];
   //               Q.q27[dM0P  ][QIN.numberOfBCnodes] = ON[dM0P  ];
   //               Q.q27[d0PP  ][QIN.numberOfBCnodes] = ON[d0PP  ];
   //               Q.q27[d0MM  ][QIN.numberOfBCnodes] = ON[d0MM  ];
   //               Q.q27[d0PM  ][QIN.numberOfBCnodes] = ON[d0PM  ];
   //               Q.q27[d0MP  ][QIN.numberOfBCnodes] = ON[d0MP  ];
   //               Q.q27[d000][QIN.numberOfBCnodes] = ON[d000];
   //               Q.q27[dPPP ][QIN.numberOfBCnodes] = ON[dPPP ];
   //               Q.q27[dMMP ][QIN.numberOfBCnodes] = ON[dMMP ];
   //               Q.q27[dPMP ][QIN.numberOfBCnodes] = ON[dPMP ];
   //               Q.q27[dMPP ][QIN.numberOfBCnodes] = ON[dMPP ];
   //               Q.q27[dPPM ][QIN.numberOfBCnodes] = ON[dPPM ];
   //               Q.q27[dMMM ][QIN.numberOfBCnodes] = ON[dMMM ];
   //               Q.q27[dPMM ][QIN.numberOfBCnodes] = ON[dPMM ];
   //               Q.q27[dMPM ][QIN.numberOfBCnodes] = ON[dMPM ];

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

   //               Q.q27[dP00   ][QIN.numberOfBCnodes] = ON[dP00   ];
   //               Q.q27[dM00   ][QIN.numberOfBCnodes] = ON[dM00   ];
   //               Q.q27[d0P0   ][QIN.numberOfBCnodes] = ON[d0P0   ];
   //               Q.q27[d0M0   ][QIN.numberOfBCnodes] = ON[d0M0   ];
   //               Q.q27[d00P   ][QIN.numberOfBCnodes] = ON[d00P   ];
   //               Q.q27[d00M   ][QIN.numberOfBCnodes] = ON[d00M   ];
   //               Q.q27[dPP0  ][QIN.numberOfBCnodes] = ON[dPP0  ];
   //               Q.q27[dMM0  ][QIN.numberOfBCnodes] = ON[dMM0  ];
   //               Q.q27[dPM0  ][QIN.numberOfBCnodes] = ON[dPM0  ];
   //               Q.q27[dMP0  ][QIN.numberOfBCnodes] = ON[dMP0  ];
   //               Q.q27[dP0P  ][QIN.numberOfBCnodes] = ON[dP0P  ];
   //               Q.q27[dM0M  ][QIN.numberOfBCnodes] = ON[dM0M  ];
   //               Q.q27[dP0M  ][QIN.numberOfBCnodes] = ON[dP0M  ];
   //               Q.q27[dM0P  ][QIN.numberOfBCnodes] = ON[dM0P  ];
   //               Q.q27[d0PP  ][QIN.numberOfBCnodes] = ON[d0PP  ];
   //               Q.q27[d0MM  ][QIN.numberOfBCnodes] = ON[d0MM  ];
   //               Q.q27[d0PM  ][QIN.numberOfBCnodes] = ON[d0PM  ];
   //               Q.q27[d0MP  ][QIN.numberOfBCnodes] = ON[d0MP  ];
   //               Q.q27[d000][QIN.numberOfBCnodes] = ON[d000];
   //               Q.q27[dPPP ][QIN.numberOfBCnodes] = ON[dPPP ];
   //               Q.q27[dMMP ][QIN.numberOfBCnodes] = ON[dMMP ];
   //               Q.q27[dPMP ][QIN.numberOfBCnodes] = ON[dPMP ];
   //               Q.q27[dMPP ][QIN.numberOfBCnodes] = ON[dMPP ];
   //               Q.q27[dPPM ][QIN.numberOfBCnodes] = ON[dPPM ];
   //               Q.q27[dMMM ][QIN.numberOfBCnodes] = ON[dMMM ];
   //               Q.q27[dPMM ][QIN.numberOfBCnodes] = ON[dPMM ];
   //               Q.q27[dMPM ][QIN.numberOfBCnodes] = ON[dMPM ];

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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
   Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
   Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
   Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
   Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
   Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
   Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
   Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
   Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
   Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
   Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
   Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
   Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
   Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
   Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
   Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
   Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
   Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
   Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
   Q.q27[d000] = &QQ[d000*sizeQ];
   Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
   Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
   Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
   Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
   Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
   Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
   Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
   Q.q27[dMPM ] = &QQ[dMPM *sizeQ];


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
               Q.q27[dP00   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dM00   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0P0   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0M0   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d00P   ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d00M   ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPM0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMP0  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dP0P  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dM0M  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dP0M  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dM0P  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d0PP  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d0MM  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0PM  ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[d0MP  ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[d000][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPPP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dMMP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dPMP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dMPP ][QIN.numberOfBCnodes] = (real)1.f;
               Q.q27[dPPM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMMM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dPMM ][QIN.numberOfBCnodes] = (real)-1.f;
               Q.q27[dMPM ][QIN.numberOfBCnodes] = (real)-1.f;
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
   //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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



////////////////////////////////////////////////////////////////////////////////
void findQPressX0(Parameter* para, int lev)
{
   // ! CAUTION ! Do not use this function!
   // As the order of the distributions was changed in July 2022, this does not work anymore.
   // https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/14
    VF_LOG_CRITICAL("findKforQPressX0() is deprecated! - see comment above for more information");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
    Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
    Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
    Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
    Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
    Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
    Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
    Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
    Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
    Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
    Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
    Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
    Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
    Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
    Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
    Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
    Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
    Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
    Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
    Q.q27[d000] = &QQ[d000*sizeQ];
    Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
    Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
    Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
    Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
    Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
    Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
    Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
    Q.q27[dMPM ] = &QQ[dMPM *sizeQ];


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
                Q.q27[dP00   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dM00   ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[d0P0   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0M0   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d00P   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d00M   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPP0  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMM0  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dPM0  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMP0  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dP0P  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dM0M  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dP0M  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dM0P  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[d0PP  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0MM  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0PM  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0MP  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d000][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPPP ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMMP ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dPMP ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMPP ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dPPM ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMMM ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dPMM ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dMPM ][QIN.numberOfBCnodes] = (real)1.f;
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
    //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
    //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
    Q.q27[dP00   ] = &QQ[dP00   *sizeQ];
    Q.q27[dM00   ] = &QQ[dM00   *sizeQ];
    Q.q27[d0P0   ] = &QQ[d0P0   *sizeQ];
    Q.q27[d0M0   ] = &QQ[d0M0   *sizeQ];
    Q.q27[d00P   ] = &QQ[d00P   *sizeQ];
    Q.q27[d00M   ] = &QQ[d00M   *sizeQ];
    Q.q27[dPP0  ] = &QQ[dPP0  *sizeQ];
    Q.q27[dMM0  ] = &QQ[dMM0  *sizeQ];
    Q.q27[dPM0  ] = &QQ[dPM0  *sizeQ];
    Q.q27[dMP0  ] = &QQ[dMP0  *sizeQ];
    Q.q27[dP0P  ] = &QQ[dP0P  *sizeQ];
    Q.q27[dM0M  ] = &QQ[dM0M  *sizeQ];
    Q.q27[dP0M  ] = &QQ[dP0M  *sizeQ];
    Q.q27[dM0P  ] = &QQ[dM0P  *sizeQ];
    Q.q27[d0PP  ] = &QQ[d0PP  *sizeQ];
    Q.q27[d0MM  ] = &QQ[d0MM  *sizeQ];
    Q.q27[d0PM  ] = &QQ[d0PM  *sizeQ];
    Q.q27[d0MP  ] = &QQ[d0MP  *sizeQ];
    Q.q27[d000] = &QQ[d000*sizeQ];
    Q.q27[dPPP ] = &QQ[dPPP *sizeQ];
    Q.q27[dMMP ] = &QQ[dMMP *sizeQ];
    Q.q27[dPMP ] = &QQ[dPMP *sizeQ];
    Q.q27[dMPP ] = &QQ[dMPP *sizeQ];
    Q.q27[dPPM ] = &QQ[dPPM *sizeQ];
    Q.q27[dMMM ] = &QQ[dMMM *sizeQ];
    Q.q27[dPMM ] = &QQ[dPMM *sizeQ];
    Q.q27[dMPM ] = &QQ[dMPM *sizeQ];


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
                Q.q27[dP00   ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dM00   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0P0   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0M0   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d00P   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d00M   ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPP0  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMM0  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPM0  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMP0  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dP0P  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dM0M  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dP0M  ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dM0P  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0PP  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0MM  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0PM  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d0MP  ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[d000][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPPP ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMMP ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPMP ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMPP ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPPM ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMMM ][QIN.numberOfBCnodes] = (real)-1.f;
                Q.q27[dPMM ][QIN.numberOfBCnodes] = (real)1.f;
                Q.q27[dMPM ][QIN.numberOfBCnodes] = (real)-1.f;
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
    //////////////  dP00   dM00   d0P0   d0M0   d00P   d00M  dPP0  dMM0  dPM0  dMP0  dP0P  dM0M  dP0M  dM0P  d0PP  d0MM  d0PM  d0MP ZERO dPPP dPPM dPMP dPMM dMPP dMPM dMMP dMMM  ////////////////////////
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
