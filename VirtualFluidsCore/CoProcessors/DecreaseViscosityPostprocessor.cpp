/*
*  DecreaseViscosityPostprocessor
*
*  Created on: 10.05.2013
*  Author: uphoff
*/

#include "DecreaseViscosityPostprocessor.h"
#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>

using namespace std;

DecreaseViscosityPostprocessor::DecreaseViscosityPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                                                               mu::Parser* nueFunc, CommunicatorPtr comm)

                                                               : Postprocessor(grid, s)
                                                               ,nueFunc(nueFunc)
                                                               ,comm(comm)
{
   if (comm->getProcessID() == comm->getRoot())
   {

   }
}
//////////////////////////////////////////////////////////////////////////
DecreaseViscosityPostprocessor::~DecreaseViscosityPostprocessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void DecreaseViscosityPostprocessor::update(double step)
{
   if(scheduler->isDue(step) )
      collectPostprocessData(step);
}
//////////////////////////////////////////////////////////////////////////
void DecreaseViscosityPostprocessor::collectPostprocessData(double step)
{

   UBLOG(logDEBUG3, "DecreaseViscosityPostprocessor::update:" << step);
   int gridRank = grid->getRank();
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   if (comm->getProcessID() == comm->getRoot())
   {

      for(int level = minInitLevel; level<=maxInitLevel;level++)
      {
         vector<Block3DPtr> blockVector;
         grid->getBlocks(level, gridRank, blockVector);
         BOOST_FOREACH(Block3DPtr block, blockVector)
         {
            LBMKernel3DPtr kernel =block->getKernel();
         }
      }

      int istep = static_cast<int>(step);
      this->timeStep       = istep;
      nueFunc->DefineVar("t" , &this->timeStep);
      double nue=nueFunc->Eval();

      for(int level = minInitLevel; level<=maxInitLevel;level++)
      {
         vector<Block3DPtr> blockVector;
         grid->getBlocks(level, gridRank, blockVector);
         BOOST_FOREACH(Block3DPtr block, blockVector)
         {
            LBMKernel3DPtr kernel =block->getKernel();
            if(kernel)      
            {
               LBMReal collFactor = LBMSystem::calcCollisionFactor(nue, block->getLevel());
               kernel->setCollisionFactor(collFactor);
            }
         }
      }

   }
}
