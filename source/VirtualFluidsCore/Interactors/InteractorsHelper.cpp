#include "InteractorsHelper.h"

#include <SetSolidBlockVisitor.h>
#include <Grid3DVisitor.h>
#include <Grid3D.h>
#include <Interactor3D.h>
#include "Block3D.h"
#include "Communicator.h"


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
void InteractorsHelper::setBC()
{
    for(Interactor3DPtr i : interactors)
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
    for(Interactor3DPtr interactor : interactors)
    {
        setBlocks(interactor, BlockType::SOLID);

        std::vector<Block3DPtr>& sb = interactor->getSolidBlockSet();
        solidBlocks.insert(solidBlocks.end(), sb.begin(), sb.end());
        interactor->removeSolidBlocks();
    }

    updateGrid();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBlocks(const Interactor3DPtr interactor, BlockType type) const
{
    SetSolidBlockVisitor v(interactor, type);
    grid->accept(v);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBcBlocks()
{
    for(const Interactor3DPtr interactor : interactors)
        setBlocks(interactor, BlockType::BC);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::updateGrid()
{
    std::vector<int> ids;

    for(const Block3DPtr block : solidBlocks)
        ids.push_back(block->getGlobalID());

    std::vector<int> rids;
    Communicator::getInstance()->allGather(ids, rids);
    grid->deleteBlocks(rids);
}

