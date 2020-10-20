#include "SetInterpolationDirsBlockVisitor.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "Grid3DSystem.h"
#include <D3Q27System.h>

SetInterpolationDirsBlockVisitor::SetInterpolationDirsBlockVisitor(std::vector<int> &dirs)
    : Block3DVisitor(0, Grid3DSystem::MAXLEVEL), dirs(dirs)
{
}
//////////////////////////////////////////////////////////////////////////
void SetInterpolationDirsBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
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
                    case NE:
                        checkFlagDir(grid, E, N, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case SW:
                        checkFlagDir(grid, W, S, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case SE:
                        checkFlagDir(grid, E, S, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case NW:
                        checkFlagDir(grid, W, N, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TE:
                        checkFlagDir(grid, E, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BW:
                        checkFlagDir(grid, W, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BE:
                        checkFlagDir(grid, E, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TW:
                        checkFlagDir(grid, W, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TN:
                        checkFlagDir(grid, N, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BS:
                        checkFlagDir(grid, S, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BN:
                        checkFlagDir(grid, N, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TS:
                        checkFlagDir(grid, S, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TNE:
                        checkFlagDir(grid, E, N, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TSW:
                        checkFlagDir(grid, W, S, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TSE:
                        checkFlagDir(grid, E, S, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case TNW:
                        checkFlagDir(grid, W, N, T, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BNE:
                        checkFlagDir(grid, E, N, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BSW:
                        checkFlagDir(grid, W, S, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BSE:
                        checkFlagDir(grid, E, S, B, flagDir, ix1, ix2, ix3, level);
                        if (!flagDir)
                            continue;
                        break;
                    case BNW:
                        checkFlagDir(grid, W, N, B, flagDir, ix1, ix2, ix3, level);
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
