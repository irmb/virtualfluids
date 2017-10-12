#include "D3Q27System.h"
namespace D3Q27System
{
   //index             0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18//falsch
   //f:              ZERO, E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
   //const int EX1[] = { 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1 };
   //const int EX2[] = { 0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1 };
   //const int EX3[] = { 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1,  1, -1, -1, -1, -1 };

   //index             0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18
   //f:                E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
   const int DX1[] = { 1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1 };
   const int DX2[] = { 0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1 };
   const int DX3[] = { 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1,  1, -1, -1, -1, -1 };

   ////index                0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18
   ////f:                   E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
   const double WEIGTH[] = { c2o27, c2o27,  c2o27,  c2o27,  c2o27,  c2o27,  c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o54, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216, c1o216 , c8o27};
  

   const int INVDIR[] = { 
                          INV_E,   
                          INV_W,  
                          INV_N,  
                          INV_S,  
                          INV_T,  
                          INV_B,  
                          INV_NE, 
                          INV_SW, 
                          INV_SE, 
                          INV_NW,
                          INV_TE, 
                          INV_BW, 
                          INV_BE, 
                          INV_TW, 
                          INV_TN, 
                          INV_BS, 
                          INV_BN, 
                          INV_TS, 
                          INV_TNE,
                          INV_TNW,
                          INV_TSE,
                          INV_TSW,
                          INV_BNE,
                          INV_BNW,
                          INV_BSE,
                          INV_BSW    };
}

    //const int FSTARTDIR = 0;
    //const int FENDDIR   = 25;   //D3Q27

    //const int STARTF = 0;
    //const int ENDF   = 26;   //D3Q27

    //const int EX1[ENDF+1];
    //const int EX2[ENDF+1];
    //const int EX3[ENDF+1];

    //const int STARTDIR = 0;
    //const int ENDDIR   = 26; //alle geometrischen richtungen

    //const int DX1[ENDDIR+1];
    //const int DX2[ENDDIR+1];
    //const int DX3[ENDDIR+1];


    //const int E    /*f1 */ = 0;
    //const int W    /*f2 */ = 1;
    //const int N    /*f3 */ = 2;
    //const int S    /*f4 */ = 3;
    //const int T    /*f5 */ = 4;
    //const int B    /*f6 */ = 5;
    //const int NE   /*f7 */ = 6;
    //const int SW   /*f8 */ = 7;
    //const int SE   /*f9 */ = 8;
    //const int NW   /*f10*/ = 9;
    //const int TE   /*f11*/ = 10;
    //const int BW   /*f12*/ = 11;
    //const int BE   /*f13*/ = 12;
    //const int TW   /*f14*/ = 13;
    //const int TN   /*f15*/ = 14;
    //const int BS   /*f16*/ = 15;
    //const int BN   /*f17*/ = 16;
    //const int TS   /*f18*/ = 17;
    //const int TNE          = 18;
    //const int TNW          = 19;
    //const int TSE          = 20;
    //const int TSW          = 21;
    //const int BNE          = 22;
    //const int BNW          = 23;
    //const int BSE          = 24;
    //const int BSW          = 25;
    //const int ZERO /*f0 */ = 26;

    //const int INV_E   = W;  
    //const int INV_W   = E;  
    //const int INV_N   = S;  
    //const int INV_S   = N;  
    //const int INV_T   = B;  
    //const int INV_B   = T;  
    //const int INV_NE  = SW; 
    //const int INV_SW  = NE; 
    //const int INV_SE  = NW; 
    //const int INV_NW  = SE; 
    //const int INV_TE  = BW; 
    //const int INV_BW  = TE; 
    //const int INV_BE  = TW; 
    //const int INV_TW  = BE; 
    //const int INV_TN  = BS; 
    //const int INV_BS  = TN; 
    //const int INV_BN  = TS; 
    //const int INV_TS  = BN; 
    //const int INV_TNE = BSW;
    //const int INV_TNW = BSE;
    //const int INV_TSE = BNW;
    //const int INV_TSW = BNE;
    //const int INV_BNE = TSW;
    //const int INV_BNW = TSE;
    //const int INV_BSE = TNW;
    //const int INV_BSW = TNE;

    //const int INVDIR[ENDDIR+1];

    //const int M_RHO     = 0;  
    //const int M_EN      = 1;  
    //const int M_EPS     = 2;  
    //const int M_JX1     = 3;  
    //const int M_QX1     = 4;  
    //const int M_JX2     = 5;  
    //const int M_QX2     = 6;  
    //const int M_JX3     = 7;  
    //const int M_QX3     = 8;  
    //const int M_3PX1X1  = 9;  
    //const int M_3PIX1X1 = 10; 
    //const int M_PWW     = 11; 
    //const int M_PIWW    = 12; 
    //const int M_PX1X2   = 13; 
    //const int M_PX2X3   = 14; 
    //const int M_PX1X3   = 15; 
    //const int M_MX1     = 16; 
    //const int M_MX2     = 17; 
    //const int M_MX3     = 18; 

    //const int STARTM = 0;
    //const int ENDM   = 18;   //D3Q27
