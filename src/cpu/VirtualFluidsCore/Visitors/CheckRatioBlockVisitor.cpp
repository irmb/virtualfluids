#include "CheckRatioBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"

CheckRatioBlockVisitor::CheckRatioBlockVisitor(int levelDepth /*shut be maxGridLevel*/, bool includeNotActiveBlocks)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), levelDepth(levelDepth), includeNotActiveBlocks(includeNotActiveBlocks),
      state(true)
{
}
//////////////////////////////////////////////////////////////////////////
void CheckRatioBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    int ix1, ix2, ix3, level;
    ix1   = block->getX1();
    ix2   = block->getX2();
    ix3   = block->getX3();
    level = block->getLevel();

    //   int nix1, nix2, nix3;
    int nlev;
    //   int neighix1, neighix2, neighix3,
    int neighlev;
    std::vector<SPtr<Block3D>> neighbors;
    grid->getAllNeighbors(ix1, ix2, ix3, level, this->levelDepth, neighbors);
    //   bool hasAdded = false;
    for (size_t i = 0; i < neighbors.size(); i++) {
        if ((neighbors[i]->isActive() || includeNotActiveBlocks) && neighbors[i]->getLevel() > level) {
            // neighix1 = neighbors[i]->getX1();
            // neighix2 = neighbors[i]->getX2();
            // neighix3 = neighbors[i]->getX3();
            neighlev = neighbors[i]->getLevel();
            // nix1 = neighix1>>1;
            // nix2 = neighix2>>1;
            // nix3 = neighix3>>1;
            nlev = neighlev - 1;

            if (nlev != level) {
                // throw UbException(UB_EXARGS, "OverlapBlockVisitor::adaptBlock - leveldifferenz passt nicht, block:
                // "+block->toString()); grid->expandBlock(ix1, ix2, ix3, level);
                state      = state && false;
                falseBlock = block;

            } else {
                state = state && true;
            }

            // UBLOG(logINFO, "OverlapBlockVisitor::state= "<<state);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
bool CheckRatioBlockVisitor::getState() { return state; }
//////////////////////////////////////////////////////////////////////////
void CheckRatioBlockVisitor::resetState() { state = true; }
//////////////////////////////////////////////////////////////////////////
std::string CheckRatioBlockVisitor::getStateString() { return falseBlock->toString(); }
//////////////////////////////////////////////////////////////////////////
