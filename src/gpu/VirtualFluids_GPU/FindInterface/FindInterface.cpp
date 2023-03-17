#include "FindInterface/FindInterface.h"

void interpolation(InterpolationCellCoarseToFine &intCF, InterpolationCellFineToCoarse &intFC, 
                   unsigned int LxCoarse, unsigned int LyCoarse, unsigned int LzCoarse, 
                   unsigned int LxFine, unsigned int LyFine, unsigned int LzFine, 
                   unsigned int dNx, unsigned int dNy, unsigned int dNz, 
                   unsigned int *kCoarse, unsigned int *kFine, bool* needInterface,
                   InterpolationCellNeighbor &offCF, InterpolationCellNeighbor &offFC)
{
   unsigned int iC,iF,jC,jF,hC,hF;
   unsigned int posCSWB, posFSWB;
   unsigned int posC;
   real xOff = (real)0.0f;
   real yOff = (real)0.0f; 
   real zOff = (real)0.0f;
   intCF.kCF    = 0;
   intFC.kFC    = 0;

   ///////////////////////////////////////////////////////////////////////////
   //Defines
   ///////////////////////////////////////////////////////////////////////////
   //Coarse to Fine
   unsigned int CF_Coarse[6], CF_NCoarse[6], CF_Fine[6], CF_NFine[6];
   unsigned int CF_xDefaultCoarse, CF_yDefaultCoarse, CF_zDefaultCoarse;
   unsigned int CF_xDefaultFine,   CF_yDefaultFine,   CF_zDefaultFine;

   CF_Coarse[ INTERFACE_E] = dNx+LxFine/2-1;
   CF_Coarse[ INTERFACE_W] = dNx;
   CF_Coarse[ INTERFACE_N] = dNy+LyFine/2-1;
   CF_Coarse[ INTERFACE_S] = dNy;
   CF_Coarse[ INTERFACE_T] = dNz+LzFine/2-1;
   CF_Coarse[ INTERFACE_B] = dNz;

   CF_NCoarse[INTERFACE_E] = dNx+LxFine/2-2;
   CF_NCoarse[INTERFACE_W] = dNx+1;
   CF_NCoarse[INTERFACE_N] = dNy+LyFine/2-2;
   CF_NCoarse[INTERFACE_S] = dNy+1;
   CF_NCoarse[INTERFACE_T] = dNz+LzFine/2-2;
   CF_NCoarse[INTERFACE_B] = dNz+1;

   CF_Fine[ INTERFACE_E]   = LxFine-2;
   CF_Fine[ INTERFACE_W]   = 0;
   CF_Fine[ INTERFACE_N]   = LyFine-2;
   CF_Fine[ INTERFACE_S]   = 0;
   CF_Fine[ INTERFACE_T]   = LzFine-2;
   CF_Fine[ INTERFACE_B]   = 0;

   CF_NFine[INTERFACE_E]   = LxFine-3;
   CF_NFine[INTERFACE_W]   = 1;
   CF_NFine[INTERFACE_N]   = LyFine-3;
   CF_NFine[INTERFACE_S]   = 1;
   CF_NFine[INTERFACE_T]   = LzFine-3;
   CF_NFine[INTERFACE_B]   = 1;

   CF_xDefaultCoarse=dNx+1;
   CF_yDefaultCoarse=dNy+1;
   CF_zDefaultCoarse=dNz+1;
   CF_xDefaultFine  =2;
   CF_yDefaultFine  =2;
   CF_zDefaultFine  =2;
   ///////////////////////////////////////////////////////////////////////////
   //Fine to Coarse
   unsigned int FC_Coarse[6], FC_NCoarse[6], FC_Fine[6], FC_NFine[6];
   unsigned int FC_xDefaultCoarse, FC_yDefaultCoarse, FC_zDefaultCoarse;
   unsigned int FC_xDefaultFine,   FC_yDefaultFine,   FC_zDefaultFine;

   FC_Coarse[ INTERFACE_E] = dNx+LxFine/2-2;
   FC_Coarse[ INTERFACE_W] = dNx+2;
   FC_Coarse[ INTERFACE_N] = dNy+LyFine/2-2;
   FC_Coarse[ INTERFACE_S] = dNy+2;
   FC_Coarse[ INTERFACE_T] = dNz+LzFine/2-2;
   FC_Coarse[ INTERFACE_B] = dNz+2;

   FC_NCoarse[INTERFACE_E] = dNx+LxFine/2-1;
   FC_NCoarse[INTERFACE_W] = dNx+1;
   FC_NCoarse[INTERFACE_N] = dNy+LyFine/2-1;
   FC_NCoarse[INTERFACE_S] = dNy+1;
   FC_NCoarse[INTERFACE_T] = dNz+LzFine/2-1;
   FC_NCoarse[INTERFACE_B] = dNz+1;

   FC_Fine[ INTERFACE_E]   = LxFine-5;
   FC_Fine[ INTERFACE_W]   = 3;
   FC_Fine[ INTERFACE_N]   = LyFine-5;
   FC_Fine[ INTERFACE_S]   = 3;
   FC_Fine[ INTERFACE_T]   = LzFine-5;
   FC_Fine[ INTERFACE_B]   = 3;

   FC_NFine[INTERFACE_E]   = LxFine-3;
   FC_NFine[INTERFACE_W]   = 1;
   FC_NFine[INTERFACE_N]   = LyFine-3;
   FC_NFine[INTERFACE_S]   = 1;
   FC_NFine[INTERFACE_T]   = LzFine-3;
   FC_NFine[INTERFACE_B]   = 1;

   FC_xDefaultCoarse=dNx+3;
   FC_yDefaultCoarse=dNy+3;
   FC_zDefaultCoarse=dNz+3;
   FC_xDefaultFine  =5;
   FC_yDefaultFine  =5;
   FC_zDefaultFine  =5;
   ///////////////////////////////////////////////////////////////////////////




   ///////////////////////////////////////////////////////////////////////////
   ///////////                                                     ///////////
   ///////////                          6                          ///////////
   ///////////                                                     ///////////
   ///////////                        FACES                        ///////////
   ///////////                                                     ///////////
   ///////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////
   //  face west                        1                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_W]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      iC = CF_Coarse[INTERFACE_W];                iF = CF_Fine[INTERFACE_W];

      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF <= LyFine-4; jC++,jF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  face north                       2                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_N]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      jC = CF_Coarse[ INTERFACE_N];               jF = CF_Fine[ INTERFACE_N];

      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {			
         for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  face east                        3                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_E]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      iC = CF_Coarse[ INTERFACE_E];               iF = CF_Fine[ INTERFACE_E];

      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF <= LyFine-4; jC++,jF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  face south                       4                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_S]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      jC = CF_Coarse[ INTERFACE_S];               jF = CF_Fine[ INTERFACE_S];

      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {			
         for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  face top                         5                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_T]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      hC = CF_Coarse[ INTERFACE_T];               hF = CF_Fine[ INTERFACE_T];

      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   =xOff;
            offFC.y[intFC.kFC]   =yOff;
            offFC.z[intFC.kFC]   =zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  face bottom                      6                                   //
   ///////////////////////////////////////////////////////////////////////////
   if (needInterface[INTERFACE_B]==true)
   {
      //////////////////////////   coarse->fine   ////////////////////////////
      hC = CF_Coarse[ INTERFACE_B];               hF = CF_Fine[ INTERFACE_B];

      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
         {
            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
            offCF.x[intCF.kCF]   = xOff;
            offCF.y[intCF.kCF]   = yOff;
            offCF.z[intCF.kCF]   = zOff;
            intCF.kCF++;
         }
      }
      //////////////////////////   fine->coarse   ////////////////////////////
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   ///////////                                                     ///////////
   ///////////                         12                          ///////////
   ///////////                                                     ///////////
   ///////////                        EDGES                        ///////////
   ///////////                                                     ///////////
   ///////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////
   //  edge east-north                  1                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || (needInterface[INTERFACE_N]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];
         iF = CF_NFine[  INTERFACE_E];
         xOff = (real)0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_E];                    
         iF = CF_Fine[   INTERFACE_E];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];
         jF = CF_NFine[  INTERFACE_N];
         yOff = (real)0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_N];                    
         jF = CF_Fine[   INTERFACE_N];
      }
      //////////////////////////////////////////////////////////////////////////
      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
      else if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge east-south                  2                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || (needInterface[INTERFACE_S]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];
         iF = CF_NFine[  INTERFACE_E];
         xOff = (real)0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_E];                    
         iF = CF_Fine[   INTERFACE_E];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];
         jF = CF_NFine[  INTERFACE_S];
         yOff = (real)-0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_S];                    
         jF = CF_Fine[   INTERFACE_S];
      }
      //////////////////////////////////////////////////////////////////////////
      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
      else if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge east-top                    3                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];
         iF = CF_NFine[  INTERFACE_E];
         xOff = (real)0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_E];                    
         iF = CF_Fine[   INTERFACE_E];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];                   
         hF = CF_NFine[  INTERFACE_T];
         zOff = (real)0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_T];                   
         hF = CF_Fine[   INTERFACE_T];
      }
      //////////////////////////////////////////////////////////////////////////
      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge east-bottom                 4                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];
         iF = CF_NFine[  INTERFACE_E];
         xOff = (real)0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_E];                    
         iF = CF_Fine[   INTERFACE_E];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];                   
         hF = CF_NFine[  INTERFACE_B];
         zOff = (real)-0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_B];                   
         hF = CF_Fine[   INTERFACE_B];
      }
      //////////////////////////////////////////////////////////////////////////
      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge west-north                  5                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || (needInterface[INTERFACE_N]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];
         iF = CF_NFine[  INTERFACE_W];
         xOff = (real)-0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_W];                    
         iF = CF_Fine[   INTERFACE_W];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];
         jF = CF_NFine[  INTERFACE_N];
         yOff = (real)0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_N];                    
         jF = CF_Fine[   INTERFACE_N];
      }
      //////////////////////////////////////////////////////////////////////////
      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
      else if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge west-south                  6                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || (needInterface[INTERFACE_S]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];
         iF = CF_NFine[  INTERFACE_W];
         xOff = (real)-0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_W];                    
         iF = CF_Fine[   INTERFACE_W];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];
         jF = CF_NFine[  INTERFACE_S];
         yOff = (real)-0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_S];                    
         jF = CF_Fine[   INTERFACE_S];
      }
      //////////////////////////////////////////////////////////////////////////
      for (hF = CF_zDefaultFine, hC = CF_zDefaultCoarse ; hF <= LzFine-4; hC++,hF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];

      for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
      else if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];

         for (hF = FC_zDefaultFine, hC = FC_zDefaultCoarse ; hF<=LzFine-7; hC++,hF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge west-top                    7                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];
         iF = CF_NFine[  INTERFACE_W];
         xOff = (real)-0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_W];                    
         iF = CF_Fine[   INTERFACE_W];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];                   
         hF = CF_NFine[  INTERFACE_T];
         zOff = (real)0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_T];                   
         hF = CF_Fine[   INTERFACE_T];
      }
      //////////////////////////////////////////////////////////////////////////
      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge west-bottom                 8                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];
         iF = CF_NFine[  INTERFACE_W];
         xOff = (real)-0.5f;
      } 
      else
      {
         iC = CF_Coarse[ INTERFACE_W];                    
         iF = CF_Fine[   INTERFACE_W];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];                   
         hF = CF_NFine[  INTERFACE_B];
         zOff = (real)-0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_B];                   
         hF = CF_Fine[   INTERFACE_B];
      }
      //////////////////////////////////////////////////////////////////////////
      for (jF = CF_yDefaultFine, jC = CF_yDefaultCoarse ; jF<=LyFine-4; jC++,jF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         for (jF = FC_yDefaultFine, jC = FC_yDefaultCoarse ; jF<=LyFine-7; jC++,jF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge north-top                   9                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_N]==true) || (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];
         jF = CF_NFine[  INTERFACE_N];
         yOff = (real)0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_N];                    
         jF = CF_Fine[   INTERFACE_N];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];                   
         hF = CF_NFine[  INTERFACE_T];
         zOff = (real)0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_T];                   
         hF = CF_Fine[   INTERFACE_T];
      }
      //////////////////////////////////////////////////////////////////////////
      for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[ INTERFACE_T]==false)
      {
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge north-bottom               10                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_N]==true) || (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];
         jF = CF_NFine[  INTERFACE_N];
         yOff = (real)0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_N];                    
         jF = CF_Fine[   INTERFACE_N];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];                   
         hF = CF_NFine[  INTERFACE_B];
         zOff = (real)-0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_B];                   
         hF = CF_Fine[   INTERFACE_B];
      }
      //////////////////////////////////////////////////////////////////////////
      for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[ INTERFACE_N]==false)
      {
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[ INTERFACE_B]==false)
      {
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge south-top                   11                                  //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_S]==true) || (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];
         jF = CF_NFine[  INTERFACE_S];
         yOff = (real)-0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_S];                    
         jF = CF_Fine[   INTERFACE_S];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];                   
         hF = CF_NFine[  INTERFACE_T];
         zOff = (real)0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_T];                   
         hF = CF_Fine[   INTERFACE_T];
      }
      //////////////////////////////////////////////////////////////////////////
      for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[ INTERFACE_T]==false)
      {
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  edge south-bottom               12                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_S]==true) || (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];
         jF = CF_NFine[  INTERFACE_S];
         yOff = (real)-0.5f;
      } 
      else
      {
         jC = CF_Coarse[ INTERFACE_S];                    
         jF = CF_Fine[   INTERFACE_S];
      }
      //////////////////////////////////////////////////////////////////////////
      if (needInterface[ INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];                   
         hF = CF_NFine[  INTERFACE_B];
         zOff = (real)-0.5f;
      } 
      else
      {
         hC = CF_Coarse[ INTERFACE_B];                   
         hF = CF_Fine[   INTERFACE_B];
      }
      //////////////////////////////////////////////////////////////////////////
      for (iF = CF_xDefaultFine, iC = CF_xDefaultCoarse ; iF<=LxFine-4; iC++,iF+=2)
      {
         posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
         intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
         offCF.x[intCF.kCF]   = xOff;
         offCF.y[intCF.kCF]   = yOff;
         offCF.z[intCF.kCF]   = zOff;
         intCF.kCF++;
      }

      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
      {			
         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[ INTERFACE_S]==false)
      {
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      } 
      else if (needInterface[ INTERFACE_B]==false)
      {
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         for (iF = FC_xDefaultFine, iC = FC_xDefaultCoarse ; iF<=LxFine-7; iC++,iF+=2)
         {			
            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
            offFC.x[intFC.kFC]   = xOff;
            offFC.y[intFC.kFC]   = yOff;
            offFC.z[intFC.kFC]   = zOff;
            intFC.kFC++;
         }
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   ///////////                                                     ///////////
   ///////////                          8                          ///////////
   ///////////                                                     ///////////
   ///////////                       CORNERS                       ///////////
   ///////////                                                     ///////////
   ///////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////
   //  corner east-north-top            1                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || 
       (needInterface[INTERFACE_N]==true) || 
       (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_E];               iF = CF_Fine[ INTERFACE_E];
      jC = CF_Coarse[ INTERFACE_N];               jF = CF_Fine[ INTERFACE_N];
      hC = CF_Coarse[ INTERFACE_T];               hF = CF_Fine[ INTERFACE_T];

      if (needInterface[INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];               iF = CF_NFine[INTERFACE_E];
         xOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];               jF = CF_NFine[INTERFACE_N];
         yOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];               hF = CF_NFine[INTERFACE_T];
         zOff = (real)0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_N]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_N]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner east-north-bottom         2                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || 
       (needInterface[INTERFACE_N]==true) || 
       (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_E];               iF = CF_Fine[ INTERFACE_E];
      jC = CF_Coarse[ INTERFACE_N];               jF = CF_Fine[ INTERFACE_N];
      hC = CF_Coarse[ INTERFACE_B];               hF = CF_Fine[ INTERFACE_B];

      if (needInterface[INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];               iF = CF_NFine[INTERFACE_E];
         xOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];               jF = CF_NFine[INTERFACE_N];
         yOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];               hF = CF_NFine[INTERFACE_B];
         zOff = (real)-0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_N]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_N]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner east-south-top            3                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || 
       (needInterface[INTERFACE_S]==true) || 
       (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_E];               iF = CF_Fine[ INTERFACE_E];
      jC = CF_Coarse[ INTERFACE_S];               jF = CF_Fine[ INTERFACE_S];
      hC = CF_Coarse[ INTERFACE_T];               hF = CF_Fine[ INTERFACE_T];

      if (needInterface[INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];               iF = CF_NFine[INTERFACE_E];
         xOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];               jF = CF_NFine[INTERFACE_S];
         yOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];               hF = CF_NFine[INTERFACE_T];
         zOff = (real)0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_S]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_S]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner east-south-bottom         4                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_E]==true) || 
       (needInterface[INTERFACE_S]==true) || 
       (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_E];               iF = CF_Fine[ INTERFACE_E];
      jC = CF_Coarse[ INTERFACE_S];               jF = CF_Fine[ INTERFACE_S];
      hC = CF_Coarse[ INTERFACE_B];               hF = CF_Fine[ INTERFACE_B];

      if (needInterface[INTERFACE_E]==false)
      {
         iC = CF_NCoarse[INTERFACE_E];               iF = CF_NFine[INTERFACE_E];
         xOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];               jF = CF_NFine[INTERFACE_S];
         yOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];               hF = CF_NFine[INTERFACE_B];
         zOff = (real)-0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_E]==false)
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_S]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_E]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_NCoarse[INTERFACE_E];               iF = FC_NFine[INTERFACE_E];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_S]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_Coarse[ INTERFACE_E];               iF = FC_Fine[ INTERFACE_E];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner west-north-top            5                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || 
       (needInterface[INTERFACE_N]==true) || 
       (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_W];               iF = CF_Fine[ INTERFACE_W];
      jC = CF_Coarse[ INTERFACE_N];               jF = CF_Fine[ INTERFACE_N];
      hC = CF_Coarse[ INTERFACE_T];               hF = CF_Fine[ INTERFACE_T];

      if (needInterface[INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];               iF = CF_NFine[INTERFACE_W];
         xOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];               jF = CF_NFine[INTERFACE_N];
         yOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];               hF = CF_NFine[INTERFACE_T];
         zOff = (real)0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_N]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_N]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner west-north-bottom         6                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || 
       (needInterface[INTERFACE_N]==true) || 
       (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_W];               iF = CF_Fine[ INTERFACE_W];
      jC = CF_Coarse[ INTERFACE_N];               jF = CF_Fine[ INTERFACE_N];
      hC = CF_Coarse[ INTERFACE_B];               hF = CF_Fine[ INTERFACE_B];

      if (needInterface[INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];               iF = CF_NFine[INTERFACE_W];
         xOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         jC = CF_NCoarse[INTERFACE_N];               jF = CF_NFine[INTERFACE_N];
         yOff = (real)0.5f;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];               hF = CF_NFine[INTERFACE_B];
         zOff = (real)-0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_N]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_N]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_N];               jF = FC_Fine[ INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_N]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_N];               jF = FC_NFine[INTERFACE_N];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner west-south-top            7                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || 
       (needInterface[INTERFACE_S]==true) || 
       (needInterface[INTERFACE_T]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_W];               iF = CF_Fine[ INTERFACE_W];
      jC = CF_Coarse[ INTERFACE_S];               jF = CF_Fine[ INTERFACE_S];
      hC = CF_Coarse[ INTERFACE_T];               hF = CF_Fine[ INTERFACE_T];

      if (needInterface[INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];               iF = CF_NFine[INTERFACE_W];
         xOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];               jF = CF_NFine[INTERFACE_S];
         yOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         hC = CF_NCoarse[INTERFACE_T];               hF = CF_NFine[INTERFACE_T];
         zOff = (real)0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_T]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_S]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_T];               hF = FC_Fine[ INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_S]==false) && (needInterface[INTERFACE_T]==false))
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_T];               hF = FC_NFine[INTERFACE_T];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

   ///////////////////////////////////////////////////////////////////////////
   //  corner west-south-bottom         8                                   //
   ///////////////////////////////////////////////////////////////////////////
   if ((needInterface[INTERFACE_W]==true) || 
       (needInterface[INTERFACE_S]==true) || 
       (needInterface[INTERFACE_B]==true))
   {
      //////////////////////////////////////////////////////////////////////////
      // coarse->fine
      //////////////////////////////////////////////////////////////////////////
      iC = CF_Coarse[ INTERFACE_W];               iF = CF_Fine[ INTERFACE_W];
      jC = CF_Coarse[ INTERFACE_S];               jF = CF_Fine[ INTERFACE_S];
      hC = CF_Coarse[ INTERFACE_B];               hF = CF_Fine[ INTERFACE_B];

      if (needInterface[INTERFACE_W]==false)
      {
         iC = CF_NCoarse[INTERFACE_W];               iF = CF_NFine[INTERFACE_W];
         xOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         jC = CF_NCoarse[INTERFACE_S];               jF = CF_NFine[INTERFACE_S];
         yOff = (real)-0.5f;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         hC = CF_NCoarse[INTERFACE_B];               hF = CF_NFine[INTERFACE_B];
         zOff = (real)-0.5f;
      }
      //////////////////////////////////////////////////////////////////////////
      posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
      intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
      intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
      offCF.x[intCF.kCF]   = xOff;
      offCF.y[intCF.kCF]   = yOff;
      offCF.z[intCF.kCF]   = zOff;
      intCF.kCF++;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      //reset
      xOff = (real)0.0f;
      yOff = (real)0.0f;
      zOff = (real)0.0f;
      //////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////
      // fine->coarse 
      //////////////////////////////////////////////////////////////////////////
      iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
      jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
      hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

      posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
      posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
      intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
      intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
      offFC.x[intFC.kFC]   = xOff;
      offFC.y[intFC.kFC]   = yOff;
      offFC.z[intFC.kFC]   = zOff;
      intFC.kFC++;

      if (needInterface[INTERFACE_W]==false)
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_S]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if (needInterface[INTERFACE_B]==false)
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_S]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_Coarse[ INTERFACE_B];               hF = FC_Fine[ INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_W]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_NCoarse[INTERFACE_W];               iF = FC_NFine[INTERFACE_W];
         jC = FC_Coarse[ INTERFACE_S];               jF = FC_Fine[ INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
      if ((needInterface[INTERFACE_S]==false) && (needInterface[INTERFACE_B]==false))
      {
         iC = FC_Coarse[ INTERFACE_W];               iF = FC_Fine[ INTERFACE_W];
         jC = FC_NCoarse[INTERFACE_S];               jF = FC_NFine[INTERFACE_S];
         hC = FC_NCoarse[INTERFACE_B];               hF = FC_NFine[INTERFACE_B];

         posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
         posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
         intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
         intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
         offFC.x[intFC.kFC]   = xOff;
         offFC.y[intFC.kFC]   = yOff;
         offFC.z[intFC.kFC]   = zOff;
         intFC.kFC++;
      }
   }

//////////////////////////////////////////////////////////////////////////
// The End
//////////////////////////////////////////////////////////////////////////
}



////////////////////////////////////////////////////////////////////////////
////old stuff
////////////////////////////////////////////////////////////////////////////
//void interpolation(InterpolationCellCF &intCF, InterpolationCellFC &intFC, 
//   unsigned int LxCoarse, unsigned int LyCoarse, unsigned int LzCoarse, 
//   unsigned int LxFine, unsigned int LyFine, unsigned int LzFine, 
//   unsigned int dNx, unsigned int dNy, unsigned int dNz, 
//   unsigned int *kCoarse, unsigned int *kFine, bool* needInterface,
//   OffsetCF &offCF, OffsetFC &offFC)
//{
//   unsigned int iC,iF,jC,jF,hC,hF;
//   unsigned int posCSWB, posFSWB;
//   unsigned int posC;
//   intCF.kCF = 0;
//   intFC.kFC = 0;
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area west                                                            //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_W]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx;                                   iF = 0;
//      jC = dNy;                                   jF = 0;
//      hC = dNz;                                   hF = 0;
//
//      for (hF = 0, hC = dNz ; hF <= LzFine-2; hC++,hF+=2)
//      {
//         for (jF = 0, jC = dNy ; jF <= LyFine-2; jC++,jF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + 2;                               iF = 3;
//      jC = dNy + 2;                               jF = 3;
//      hC = dNz + 2;                               hF = 3;
//
//      for (hF = 3, hC = dNz + 2 ; hF<=LzFine-5; hC++,hF+=2)
//      {			
//         for (jF = 3, jC = dNy + 2 ; jF<=LyFine-5; jC++,jF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area north                                                           //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_N]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx;                                   iF = 0;
//      jC = dNy + LyFine/2 - 1;                    jF = LyFine - 2;
//      hC = dNz;                                   hF = 0;
//
//      for (hF = 0, hC = dNz ; hF <= LzFine-2; hC++,hF+=2)
//      {			
//         for (iF = 0 , iC = dNx ; iF<=LxFine-2; iC++,iF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + 2;                               iF = 3;
//      jC = dNy + LyFine/2 - 2;                    jF = LyFine - 5;
//      hC = dNz + 2;                               hF = 3;
//
//      for (hF = 3, hC = dNz + 2 ; hF<=LzFine-5; hC++,hF+=2)
//      {			
//         for (iF = 3, iC = dNx + 2 ; iF<=LxFine-5; iC++,iF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area east                                                            //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_E]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx+LxFine/2-1;                        iF = LxFine-2;
//      jC = dNy;                                   jF = 0;
//      hC = dNz;                                   hF = 0;
//
//      for (hF = 0, hC = dNz ; hF<=LzFine-2; hC++,hF+=2)
//      {			
//         for (jF = 0, jC = dNy ; jF<=LyFine-2; jC++,jF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + LxFine/2 - 2;   iF = LxFine - 5;
//      jC = dNy + 2;              jF = 3;
//      hC = dNz + 2;              hF = 3;
//
//      for (hF = 3, hC = dNz + 2 ; hF<=LzFine-5; hC++,hF+=2)
//      {			
//         for (jF = 3, jC = dNy + 2 ; jF<=LyFine-5; jC++,jF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area south                                                           //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_S]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx;                                   iF = 0;
//      jC = dNy;                                   jF = 0;
//      hC = dNz;                                   hF = 0;
//
//      for (hF = 0, hC = dNz ; hF <= LzFine-2; hC++,hF+=2)
//      {			
//         for (iF = 0, iC = dNx ; iF<=LxFine-2; iC++,iF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + 2;                               iF = 3;
//      jC = dNy + 2;                               jF = 3;
//      hC = dNz + 2;                               hF = 3;
//
//      for (hF = 3, hC = dNz + 2 ; hF<=LzFine-5; hC++,hF+=2)
//      {			
//         for (iF = 3, iC = dNx + 2 ; iF<=LxFine-5; iC++,iF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area top                                                             //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_T]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx;                                   iF = 0;
//      jC = dNy;                                   jF = 0;
//      hC = dNz;                                   hF = 0;
//
//      for (jF = 0, jC = dNy ; jF<=LyFine-2; jC++,jF+=2)
//      {
//         for (iF = 0, iC = dNx ; iF<=LxFine-2; iC++,iF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + 2;                               iF = 3;
//      jC = dNy + 2;                               jF = 3;
//      hC = dNz + 2;                               hF = 3;
//
//      for (jF = 3, jC = dNy + 2 ; jF<=LyFine-5; jC++,jF+=2)
//      {			
//         for (iF = 3, iC = dNx + 2 ; iF<=LxFine-5; iC++,iF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//
//   ///////////////////////////////////////////////////////////////////////////
//   //  area bottom                                                          //
//   ///////////////////////////////////////////////////////////////////////////
//   if (needInterface[INTERFACE_B]==true)
//   {
//      //////////////////////////   coarse->fine   ////////////////////////////
//      iC = dNx;                                   iF = 0;
//      jC = dNy;                                   jF = 0;
//      hC = dNz + LzFine/2 -1;                     hF = LzFine - 2;
//
//      for (jF = 0, jC = dNy ; jF<=LyFine-2; jC++,jF+=2)
//      {			
//         for (iF = 0, iC = dNx ; iF<=LxFine-2; iC++,iF+=2)
//         {
//            posCSWB=vectorPosition(iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine  , LyFine);
//            intCF.ICellCFC[intCF.kCF] = kCoarse[posCSWB];
//            intCF.ICellCFF[intCF.kCF] = kFine[posFSWB];
//            offCF.xOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.yOffCF[intCF.kCF]   = (real)0.0f;
//            offCF.zOffCF[intCF.kCF]   = (real)0.0f;
//            intCF.kCF++;
//         }
//      }
//      //////////////////////////   fine->coarse   ////////////////////////////
//      iC = dNx + 2;                               iF = 3; 
//      jC = dNy + 2;                               jF = 3;
//      hC = dNz + LzFine/2 - 2;                    hF = LzFine - 5;
//
//      for (jF = 3, jC = dNy + 2 ; jF<=LyFine-5; jC++,jF+=2)
//      {			
//         for (iF = 3, iC = dNx + 2 ; iF<=LxFine-5; iC++,iF+=2)
//         {			
//            posC=vectorPosition(   iC, jC, hC, LxCoarse, LyCoarse);
//            posFSWB=vectorPosition(iF, jF, hF, LxFine,   LyFine);
//            intFC.ICellFCC[intFC.kFC] = kCoarse[posC];
//            intFC.ICellFCF[intFC.kFC] = kFine[posFSWB];
//            offFC.xOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.yOffFC[intFC.kFC]   = (real)0.0f;
//            offFC.zOffFC[intFC.kFC]   = (real)0.0f;
//            intFC.kFC++;
//         }
//      }
//   }
//}
//
//
//
