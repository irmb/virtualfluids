#include "ChangeBoundaryDensityBlockVisitor.h"
#include "LBMKernelETD3Q27.h"
#include "Grid3DSystem.h"
#include "D3Q27BoundaryCondition.h"

ChangeBoundaryDensityBlockVisitor::ChangeBoundaryDensityBlockVisitor(float oldBoundaryDensity, float newBoundaryDensity) :
Block3DVisitor(0, Grid3DSystem::MAXLEVEL), 
oldBoundaryDensity(oldBoundaryDensity), 
newBoundaryDensity(newBoundaryDensity)
{

}
//////////////////////////////////////////////////////////////////////////
ChangeBoundaryDensityBlockVisitor::~ChangeBoundaryDensityBlockVisitor()
{

}
//////////////////////////////////////////////////////////////////////////
void ChangeBoundaryDensityBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if (block->getRank() == grid->getRank())
   {
      LBMKernel3DPtr kernel = block->getKernel();
      BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();

      int minX1 = 0;
      int minX2 = 0;
      int minX3 = 0;
      int maxX1 = (int)bcArray.getNX1();
      int maxX2 = (int)bcArray.getNX2();
      int maxX3 = (int)bcArray.getNX3();

      for (int x3 = minX3; x3 < maxX3; x3++)
      {
         for (int x2 = minX2; x2 < maxX2; x2++)
         {
            for (int x1 = minX1; x1 < maxX1; x1++)
            {
               if (!bcArray.isSolid(x1, x2, x3) && !bcArray.isUndefined(x1, x2, x3))
               {
                  bcPtr = bcArray.getBC(x1, x2, x3);
                  if (bcPtr)
                  {
                     if (bcPtr->hasDensityBoundary())
                     {
                        float bcDensity = bcPtr->getBoundaryDensity();
                        if (bcDensity == oldBoundaryDensity)
                        {
                           bcPtr->setBoundaryDensity(newBoundaryDensity);
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
