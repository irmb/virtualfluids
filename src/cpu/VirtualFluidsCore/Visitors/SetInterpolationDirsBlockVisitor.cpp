#include "SetInterpolationDirsBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include <D3Q27System.h>

SetInterpolationDirsBlockVisitor::SetInterpolationDirsBlockVisitor(std::vector<int> &dirs)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), dirs(dirs)
{
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationDirsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    int ix1, ix2, ix3, level;
    ix1   = block->getX1();
    ix2   = block->getX2();
    ix3   = block->getX3();
    level = block->getLevel();
    using namespace D3Q27System;
    if (level == 0)
        return;

    SPtr<Block3D> parentblock = grid->getSuperBlock(ix1, ix2, ix3, level);
    if (!parentblock)
        return;

    for (int dir : dirs) {
        SPtr<Block3D> nblock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);
        if (!nblock) {
            SPtr<Block3D> p_nblock = grid->getNeighborBlock(dir, parentblock);

            if (p_nblock) {
                bool flagDir;
                switch (dir) {
                    case DIR_PP0:
                        checkFlagDir(grid, dP00, DIR_0P0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MM0:
                        checkFlagDir(grid, dM00, DIR_0M0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_PM0:
                        checkFlagDir(grid, dP00, DIR_0M0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MP0:
                        checkFlagDir(grid, dM00, DIR_0P0, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_P0P:
                        checkFlagDir(grid, dP00, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_M0M:
                        checkFlagDir(grid, dM00, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_P0M:
                        checkFlagDir(grid, dP00, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_M0P:
                        checkFlagDir(grid, dM00, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_0PP:
                        checkFlagDir(grid, DIR_0P0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_0MM:
                        checkFlagDir(grid, DIR_0M0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_0PM:
                        checkFlagDir(grid, DIR_0P0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_0MP:
                        checkFlagDir(grid, DIR_0M0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_PPP:
                        checkFlagDir(grid, dP00, DIR_0P0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MMP:
                        checkFlagDir(grid, dM00, DIR_0M0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_PMP:
                        checkFlagDir(grid, dP00, DIR_0M0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MPP:
                        checkFlagDir(grid, dM00, DIR_0P0, DIR_00P, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_PPM:
                        checkFlagDir(grid, dP00, DIR_0P0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MMM:
                        checkFlagDir(grid, dM00, DIR_0M0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_PMM:
                        checkFlagDir(grid, dP00, DIR_0M0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case DIR_MPM:
                        checkFlagDir(grid, dM00, DIR_0P0, DIR_00M, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                }

                block->setInterpolationFlagFC(dir);
                parentblock->setInterpolationFlagCF(dir);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, bool &flagDirection, int
// ix1, int ix2, int ix3, int level)
//{
//   SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
//   SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);
//   if (!((block1 && block2)  ||  (!block1 && !block2)))
//      flagDirection = false;
//   else
//      flagDirection = true;
//}

void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, bool &flagDirection, int ix1,
                                                    int ix2, int ix3, int level)
{
    SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
    SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);

    SPtr<Block3D> pblock  = grid->getSuperBlock(ix1, ix2, ix3, level);
    SPtr<Block3D> pblock1 = grid->getNeighborBlock(dir1, pblock);
    SPtr<Block3D> pblock2 = grid->getNeighborBlock(dir2, pblock);

    if (!((block1 && block2) || (!block1 && !block2)) || !((pblock1 && pblock2) || (!pblock1 && !pblock2)))
        flagDirection = false;
    else
        flagDirection = true;
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationDirsBlockVisitor::checkFlagDir(SPtr<Grid3D> grid, int dir1, int dir2, int dir3,
                                                    bool &flagDirection, int ix1, int ix2, int ix3, int level)
{
    SPtr<Block3D> block1 = grid->getNeighborBlock(dir1, ix1, ix2, ix3, level);
    SPtr<Block3D> block2 = grid->getNeighborBlock(dir2, ix1, ix2, ix3, level);
    SPtr<Block3D> block3 = grid->getNeighborBlock(dir3, ix1, ix2, ix3, level);
    if (!((block1 && block2 && block3) || (!block1 && !block2 && !block3)))
        flagDirection = false;
    else
        flagDirection = true;
}
