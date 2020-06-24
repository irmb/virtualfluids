#include <Grid3DSystem.h>

namespace Grid3DSystem
{
   const int INVDIR[] = { INV_E  ,   
                          INV_W  ,  
                          INV_N  ,  
                          INV_S  ,  
                          INV_T  ,  
                          INV_B  ,  
                          INV_NE , 
                          INV_NW , 
                          INV_SE , 
                          INV_SW ,
                          INV_TE , 
                          INV_TW , 
                          INV_BE , 
                          INV_BW , 
                          INV_TN , 
                          INV_TS , 
                          INV_BN , 
                          INV_BS , 
                          INV_TNE,
                          INV_TNW,
                          INV_TSE,
                          INV_TSW,
                          INV_BNE,
                          INV_BNW,
                          INV_BSE,
                          INV_BSW    };

   //index             0   1   2   3   4   5  6   7   8    9  10  11  12  13  14  15  16  17  18
   //direction:        E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN, BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW
   const int EX1[] = { 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
   const int EX2[] = { 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1 };
   const int EX3[] = { 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1 };
}

//////////////////////////////////////////////////////////////////////////
const int& Grid3DSystem::getInvertDirection(const int& direction)
{  
#ifdef _DEBUG
   if(direction<STARTDIR || direction>ENDDIR) 
      throw UbException(UB_EXARGS,"unknown direction");
#endif
   return INVDIR[direction];
}

