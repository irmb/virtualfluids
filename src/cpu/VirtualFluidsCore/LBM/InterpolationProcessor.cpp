#include "InterpolationProcessor.h"


//////////////////////////////////////////////////////////////////////////
InterpolationProcessor::InterpolationProcessor()
= default;
//////////////////////////////////////////////////////////////////////////
InterpolationProcessor::~InterpolationProcessor()
= default;
//////////////////////////////////////////////////////////////////////////
void InterpolationProcessor::readICell(SPtr<DistributionArray3D> f, D3Q27ICell& icell, int x1, int x2, int x3) 
{
   f->getDistribution(icell.BSW, x1, x2, x3);
   f->getDistribution(icell.BSE, x1+1, x2, x3);
   f->getDistribution(icell.BNW, x1, x2+1, x3);
   f->getDistribution(icell.BNE, x1+1, x2+1, x3);
   f->getDistribution(icell.TSW, x1, x2, x3+1);
   f->getDistribution(icell.TSE, x1+1, x2, x3+1);
   f->getDistribution(icell.TNW, x1, x2+1, x3+1);
   f->getDistribution(icell.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationProcessor::writeICell(SPtr<DistributionArray3D> f, const D3Q27ICell& icell, int x1, int x2, int x3)
{
   f->setDistribution(icell.BSW, x1, x2, x3);
   f->setDistribution(icell.BSE, x1+1, x2, x3);
   f->setDistribution(icell.BNW, x1, x2+1, x3);
   f->setDistribution(icell.BNE, x1+1, x2+1, x3);
   f->setDistribution(icell.TSW, x1, x2, x3+1);
   f->setDistribution(icell.TSE, x1+1, x2, x3+1);
   f->setDistribution(icell.TNW, x1, x2+1, x3+1);
   f->setDistribution(icell.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationProcessor::writeICellInv(SPtr<DistributionArray3D> f, const D3Q27ICell& icell, int x1, int x2, int x3) 
{
   f->setDistributionInv(icell.BSW, x1, x2, x3);
   f->setDistributionInv(icell.BSE, x1+1, x2, x3);
   f->setDistributionInv(icell.BNW, x1, x2+1, x3);
   f->setDistributionInv(icell.BNE, x1+1, x2+1, x3);
   f->setDistributionInv(icell.TSW, x1, x2, x3+1);
   f->setDistributionInv(icell.TSE, x1+1, x2, x3+1);
   f->setDistributionInv(icell.TNW, x1, x2+1, x3+1);
   f->setDistributionInv(icell.TNE, x1+1, x2+1, x3+1);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationProcessor::writeINode(SPtr<DistributionArray3D> f, const LBMReal* const inode, int x1, int x2, int x3)
{
   f->setDistribution(inode, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
void InterpolationProcessor::writeINodeInv(SPtr<DistributionArray3D> f, const LBMReal* const inode, int x1, int x2, int x3) 
{
   f->setDistributionInv(inode, x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
bool InterpolationProcessor::iCellHasSolid(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3) 
{
   for (int ix3 = x3; ix3 <= x3 + 1; ix3++)
      for(int ix2 = x2; ix2 <= x2 + 1; ix2++)
         for(int ix1 = x1; ix1 <= x1 + 1; ix1++)
         {
            if(bcArray->isSolid(ix1, ix2, ix3))
               return true;
         }
   return false;  
}
//////////////////////////////////////////////////////////////////////////
bool InterpolationProcessor::findNeighborICell(const SPtr<BCArray3D> bcArray, SPtr<DistributionArray3D> f, 
                                                    D3Q27ICell& icell, int maxX1, int maxX2, int maxX3, 
                                                    int x1, int x2, int x3, LBMReal& xoff, LBMReal& yoff, LBMReal& zoff) 
{
   m_maxX1 = maxX1;
   m_maxX2 = maxX2;
   m_maxX3 = maxX3;

   //GoWest
   if(inRange(x1-1,x2,x3) && !iCellHasSolid(bcArray, x1-1,x2,x3))
   {
      readICell(f,icell,x1-1,x2,x3);
      xoff = 1;
      yoff = 0;
      zoff = 0;
   }
   //GoEast
   else if(inRange(x1+2,x2,x3) && !iCellHasSolid(bcArray, x1+1,x2,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2,x3);
      xoff = -1;
      yoff = 0;
      zoff = 0;
   }
   //GoSouth
   else if(inRange(x1,x2-1,x3) && !iCellHasSolid(bcArray, x1,x2-1,x3)) 
   {
      readICell(f,icell,x1,x2-1,x3);
      xoff = 0;
      yoff = 1;
      zoff = 0;
   }
   //GoNorth
   else if(inRange(x1,x2+2,x3) && !iCellHasSolid(bcArray, x1,x2+1,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1,x2+1,x3);
      xoff = 0;
      yoff = -1;
      zoff = 0;
   }
   //GoBottom
   else if(inRange(x1,x2,x3-1) && !iCellHasSolid(bcArray, x1,x2,x3-1)) 
   {
      readICell(f,icell,x1,x2,x3-1);
      xoff = 0;
      yoff = 0;
      zoff = 1;
   }
   //GoTop
   else if(inRange(x1,x2,x3+2) && !iCellHasSolid(bcArray, x1,x2,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1,x2,x3+1);
      xoff = 0;
      yoff = 0;
      zoff = -1;
   }
   //GoNW
   else if(inRange(x1-1,x2+2,x3) && !iCellHasSolid(bcArray, x1-1,x2+1,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2+1,x3);
      xoff = 1;
      yoff = -1;
      zoff = 0;
   }
   //GoNE
   else if(inRange(x1+2,x2+2,x3) && !iCellHasSolid(bcArray, x1+1,x2+1,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2+1,x3);
      xoff = -1;
      yoff = -1;
      zoff = 0;
   }
   //GoSW
   else if(inRange(x1-1,x2-1,x3) && !iCellHasSolid(bcArray, x1-1,x2-1,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2-1,x3);
      xoff = 1;
      yoff = 1;
      zoff = 0;
   }
   //GoSE
   else if(inRange(x1+2,x2-1,x3) && !iCellHasSolid(bcArray, x1+1,x2-1,x3)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2-1,x3);
      xoff = -1;
      yoff = 1;
      zoff = 0;
   }
   //GoBW
   else if(inRange(x1-1,x2,x3-1) && !iCellHasSolid(bcArray, x1-1,x2,x3-1))
   {
      readICell(f,icell,x1-1,x2,x3-1);
      xoff = 1;
      yoff = 0;
      zoff = 1;
   }
   //GoBE
   else if(inRange(x1+2,x2,x3-1) && !iCellHasSolid(bcArray, x1+1,x2,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2,x3-1);
      xoff = -1;
      yoff = 0;
      zoff = 1;
   }
   //GoBS
   else if(inRange(x1,x2-1,x3-1) && !iCellHasSolid(bcArray, x1,x2-1,x3-1)) 
   {
      readICell(f,icell,x1,x2-1,x3-1);
      xoff = 0;
      yoff = 1;
      zoff = 1;
   }
   //GoBN
   else if(inRange(x1,x2+2,x3-1) && !iCellHasSolid(bcArray, x1,x2+1,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1,x2+1,x3-1);
      xoff = 0;
      yoff = -1;
      zoff = 1;
   }
   //GoTW
   else if(inRange(x1-1,x2,x3+2) && !iCellHasSolid(bcArray, x1-1,x2,x3+1))
   {
      readICell(f,icell,x1-1,x2,x3+1);
      xoff = 1;
      yoff = 0;
      zoff = -1;
   }
   //GoTE
   else if(inRange(x1+2,x2,x3+2) && !iCellHasSolid(bcArray, x1+1,x2,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2,x3+1);
      xoff = -1;
      yoff = 0;
      zoff = -1;
   }
   //GoTS
   else if(inRange(x1,x2-1,x3+2) && !iCellHasSolid(bcArray, x1,x2-1,x3+1)) 
   {
      readICell(f,icell,x1,x2-1,x3+1);
      xoff = 0;
      yoff = 1;
      zoff = -1;
   }
   //GoTN
   else if(inRange(x1,x2+2,x3+2) && !iCellHasSolid(bcArray, x1,x2+1,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1,x2+1,x3+1);
      xoff = 0;
      yoff = -1;
      zoff = -1;
   }
   //GoTNW
   else if(inRange(x1-1,x2+2,x3+2) && !iCellHasSolid(bcArray, x1-1,x2+1,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2+1,x3+1);
      xoff = 1;
      yoff = -1;
      zoff = -1;
   }
   //GoTNE
   else if(inRange(x1+2,x2+2,x3+2) && !iCellHasSolid(bcArray, x1+1,x2+1,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2+1,x3+1);
      xoff = -1;
      yoff = -1;
      zoff = -1;
   }
   //GoTSE
   else if(inRange(x1+2,x2-1,x3+2) && !iCellHasSolid(bcArray, x1+1,x2-1,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2-1,x3+1);
      xoff = -1;
      yoff =  1;
      zoff = -1;
   }
   //GoTSW
   else if(inRange(x1-1,x2-1,x3+2) && !iCellHasSolid(bcArray, x1-1,x2-1,x3+1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2-1,x3+1);
      xoff =  1;
      yoff =  1;
      zoff = -1;
   }
   //GoBNW
   else if(inRange(x1-1,x2+2,x3-1) && !iCellHasSolid(bcArray, x1-1,x2+1,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2+1,x3-1);
      xoff =  1;
      yoff = -1;
      zoff =  1;
   }
   //GoBNE
   else if(inRange(x1+2,x2+2,x3-1) && !iCellHasSolid(bcArray, x1+1,x2+1,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2+1,x3-1);
      xoff = -1;
      yoff = -1;
      zoff =  1;
   }
   //GoBSE
   else if(inRange(x1+2,x2-1,x3-1) && !iCellHasSolid(bcArray, x1+1,x2-1,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1+1,x2-1,x3-1);
      xoff = -1;
      yoff =  1;
      zoff =  1;
   }
   //GoBSW
   else if(inRange(x1-1,x2-1,x3-1) && !iCellHasSolid(bcArray, x1-1,x2-1,x3-1)) // ist �bern�chster Knoten auch im Gebiet (Grundknoten bei 0,0,0
   {
      readICell(f,icell,x1-1,x2-1,x3-1);
      xoff =  1;
      yoff =  1;
      zoff =  1;
   }
   //default
   else
   {
      //std::string err = "For x1="+StringUtil::toString(x1)+", x2=" + StringUtil::toString(x2)+", x3=" + StringUtil::toString(x3)+
      //                  " interpolation is not implemented for other direction"+
      //                  " by using in: "+(std::string)typeid(*this).name()+ 
      //                  " or maybe you have a solid on the block boundary";
      //UB_THROW(UbException(UB_EXARGS, err));
      return false;
   }
   return true;
}
//////////////////////////////////////////////////////////////////////////
int InterpolationProcessor::iCellHowManySolids( const SPtr<BCArray3D> bcArray, int x1, int x2, int x3 )
{
   int count = 0;
   for (int ix3 = x3; ix3 <= x3 + 1; ix3++)
      for(int ix2 = x2; ix2 <= x2 + 1; ix2++)
         for(int ix1 = x1; ix1 <= x1 + 1; ix1++)
         {
            if(bcArray->isSolid(ix1, ix2, ix3))
               count++;
         }
   return count;  
}
