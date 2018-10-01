#include "InteractorsHelper.h"

#include <Grid3DVisitor.h>
#include <Grid3D.h>
#include <Interactor3D.h>
#include "Block3D.h"
#include "Communicator.h"
#include "SetSolidBlocksBlockVisitor.h"
#include "SetBcBlocksBlockVisitor.h"


InteractorsHelper::InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor, bool deleteBlocks) :
                                     grid(grid), visitor(visitor), deleteBlocks(deleteBlocks)
{

}
//////////////////////////////////////////////////////////////////////////
InteractorsHelper::~InteractorsHelper()
{

}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::addInteractor( SPtr<Interactor3D> interactor )
{
   interactors.push_back(interactor);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBC()
{
    for(SPtr<Interactor3D> i : interactors)
        i->initInteractor();
}

void InteractorsHelper::sendDomainDecompositionVisitor() const
{
    grid->accept( visitor );
}

//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::selectBlocks()
{
   sendDomainDecompositionVisitor();
   deleteSolidBlocks();

   sendDomainDecompositionVisitor();
   setBcBlocks();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::deleteSolidBlocks()
{
    for(SPtr<Interactor3D> interactor : interactors)
    {
        SetSolidBlocksBlockVisitor v(interactor);
        grid->accept(v);
        if (deleteBlocks)
        {
           std::vector<SPtr<Block3D>>& sb = interactor->getSolidBlockSet();
           solidBlocks.insert(solidBlocks.end(), sb.begin(), sb.end());
           interactor->removeSolidBlocks();
        }
    }

   if (deleteBlocks) updateGrid();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBcBlocks()
{
    for(const SPtr<Interactor3D> interactor : interactors)
    {
       SetBcBlocksBlockVisitor v(interactor);
       grid->accept(v);
    }
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::updateGrid()
{
    std::vector<int> ids;

    for(const SPtr<Block3D> block : solidBlocks)
        ids.push_back(block->getGlobalID());

    std::vector<int> rids;
    Communicator::getInstance()->allGather(ids, rids);
    grid->deleteBlocks(rids);
}

