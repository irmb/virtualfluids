#include "InteractorsHelper.h"

#include <SetSolidBlockVisitor.h>
#include <Grid3DVisitor.h>
#include <Grid3D.h>
#include <Interactor3D.h>
#include "Block3D.h"
#include "Communicator.h"


InteractorsHelper::InteractorsHelper(SPtr<Grid3D> grid, SPtr<Grid3DVisitor> visitor) :
                                     grid(grid), visitor(visitor)
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
        setBlocks(interactor, BlockType::SOLID);

        std::vector<SPtr<Block3D>>& sb = interactor->getSolidBlockSet();
        solidBlocks.insert(solidBlocks.end(), sb.begin(), sb.end());
        interactor->removeSolidBlocks();
    }

    updateGrid();
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBlocks(const SPtr<Interactor3D> interactor, BlockType type) const
{
    SetSolidBlockVisitor v(interactor, type);
    grid->accept(v);
}
//////////////////////////////////////////////////////////////////////////
void InteractorsHelper::setBcBlocks()
{
    for(const SPtr<Interactor3D> interactor : interactors)
        setBlocks(interactor, BlockType::BC);
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

