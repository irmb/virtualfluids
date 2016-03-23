#include "InteractorsHelper.h"
#include <SetSolidOrTransBlockVisitor.h>
#include <Grid3DVisitor.h>
#include <boost/foreach.hpp>


InteractorsHelper::InteractorsHelper(Grid3DPtr grid, Grid3DVisitorPtr visitor) :
                                     grid(grid), visitor(visitor)
{

}
//////////////////////////////////////////////////////////////////////////
InteractorsHelper::~InteractorsHelper()
{

}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::addInteractor( Interactor3DPtr interactor )
{
   interactors.push_back(interactor);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::deleteSolidBlocks()
{
   BOOST_FOREACH(Interactor3DPtr i, interactors)
   {
      //UBLOG(logINFO,"rank="<<grid->getRank()<<", SetSolidOrTransBlockVisitor::start");
      SetSolidOrTransBlockVisitor v(i, SetSolidOrTransBlockVisitor::SOLID);
      grid->accept(v);
      //UBLOG(logINFO,"rank="<<grid->getRank()<<", SetSolidOrTransBlockVisitor::end");
      std::vector<Block3DPtr>& sb = i->getSolidBlockSet();
      solidBlocks.insert(solidBlocks.end(), sb.begin(), sb.end());
      i->removeSolidBlocks();
   }
   
   //UBLOG(logINFO,"rank="<<grid->getRank()<<", solidBlocks.size = " <<solidBlocks.size());
   
   updateGrid();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::updateGrid()
{
   std::vector<int> ids;

   BOOST_FOREACH(Block3DPtr block, solidBlocks)
   {
      ids.push_back(block->getGlobalID());
   }
   
   //UBLOG(logINFO,"rank="<<grid->getRank()<<", ids.size = " <<ids.size());

   std::vector<int> rids;
   Communicator::getInstance()->allGather(ids, rids);
   grid->deleteBlocks(rids);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setTransBlocks()
{
   BOOST_FOREACH(Interactor3DPtr i, interactors)
   {
      SetSolidOrTransBlockVisitor v(i, SetSolidOrTransBlockVisitor::TRANS);
      grid->accept(v);
   }
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::selectBlocks()
{
   //domain decomposition visitor
   grid->accept( visitor );
   //delete solid blocks
   deleteSolidBlocks();
   //domain decomposition visitor
   grid->accept( visitor );
   //set trans blocks
   setTransBlocks();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBC()
{
   BOOST_FOREACH(Interactor3DPtr i, interactors)
   {
      i->initInteractor();
   }
}
