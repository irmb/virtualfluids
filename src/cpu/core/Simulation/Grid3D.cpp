//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Grid.cpp
//! \ingroup Simulation
//! \author Konstantin Kutscher
//=======================================================================================

#include "Grid3D.h"

#include <set>

#include <basics/writer/WbWriterVtkXmlASCII.h>
#include <geometry3d/CoordinateTransformation3D.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Grid3DVisitor.h"
#include "Interactor3D.h"
#include "D3Q27System.h"
#include <Block3D.h>
#include <parallel/Communicator.h>
#include "UbMath.h"

using namespace std;

Grid3D::Grid3D() { levelSet.resize(D3Q27System::MAXLEVEL + 1); }
//////////////////////////////////////////////////////////////////////////
Grid3D::Grid3D(std::shared_ptr<vf::parallel::Communicator> comm)

{
    levelSet.resize(D3Q27System::MAXLEVEL + 1);
    bundle = comm->getBundleID();
    rank = comm->getProcessID();
}
//////////////////////////////////////////////////////////////////////////
Grid3D::Grid3D(std::shared_ptr<vf::parallel::Communicator> comm, int blockNx1, int blockNx2, int blockNx3, int gridNx1, int gridNx2, int gridNx3)
    :

      blockNx1(blockNx1), blockNx2(blockNx2), blockNx3(blockNx2), nx1(gridNx1), nx2(gridNx2), nx3(gridNx3)
{
    using namespace vf::basics::constant;

    levelSet.resize(D3Q27System::MAXLEVEL + 1);
    bundle = comm->getBundleID();
    rank  = comm->getProcessID();
    trafo = std::make_shared<CoordinateTransformation3D>(c0o1, c0o1, c0o1, (real)blockNx1, (real)blockNx2,
                                                         (real)blockNx3);
    UbTupleInt3 minInd(0, 0, 0);
    UbTupleInt3 maxInd(gridNx1, gridNx2, gridNx3);
    this->fillExtentWithBlocks(minInd, maxInd);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::addInteractor(SPtr<Interactor3D> interactor) { interactors.push_back(interactor); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::addAndInitInteractor(SPtr<Interactor3D> interactor, real timestep)
{
    interactors.push_back(interactor);
    interactor->initInteractor(timestep);
}
//////////////////////////////////////////////////////////////////////////
Grid3D::Interactor3DSet Grid3D::getInteractors() { return interactors; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::accept(Block3DVisitor &blockVisitor)
{
    int startLevel = blockVisitor.getStartLevel();
    int stopLevel  = blockVisitor.getStopLevel();

    if (startLevel < 0 || stopLevel < 0 || startLevel > D3Q27System::MAXLEVEL || stopLevel > D3Q27System::MAXLEVEL)
        throw UbException(UB_EXARGS, "not valid level!");

    bool dir = startLevel < stopLevel;
    if (dir)
        stopLevel += 1;
    else
        stopLevel -= 1;

    //#pragma omp parallel
    //   {
    //      for (int l = startLevel; l!=stopLevel;)
    //      {
    //         std::vector<SPtr<Block3D>> blockVector;
    //         getBlocks(l, blockVector);
    //         int sizeb = (int)blockVector.size();
    //
    //#pragma omp for
    //         for (int i = 0; i < sizeb; i++)
    //         {
    //            blockVisitor.visit(shared_from_this(), blockVector[i]);
    //         }
    //         if (dir)  l++;
    //         else     l--;
    //      }
    //   }
    for (int l = startLevel; l != stopLevel;) {
        std::vector<SPtr<Block3D>> blockVector;
        getBlocks(l, blockVector);
        for (SPtr<Block3D> b : blockVector) {
            blockVisitor.visit(shared_from_this(), b);
        }
        if (dir)
            l++;
        else
            l--;
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::accept(Grid3DVisitor &gridVisitor) { gridVisitor.visit(shared_from_this()); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::accept(SPtr<Grid3DVisitor> gridVisitor) { gridVisitor->visit(shared_from_this()); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::addBlock(SPtr<Block3D> block)
{
    if (block) {
        this->blockIdMap.insert(std::make_pair(block->getGlobalID(), block));
        int level = block->getLevel();
        this->levelSet[level].insert(std::make_pair(Block3DKey(block->getX1(), block->getX2(), block->getX3()), block));
    }
}
//////////////////////////////////////////////////////////////////////////
bool Grid3D::deleteBlock(SPtr<Block3D> block)
{
    return this->deleteBlock(block->getX1(), block->getX2(), block->getX3(), block->getLevel());
}
//////////////////////////////////////////////////////////////////////////
bool Grid3D::deleteBlock(int ix1, int ix2, int ix3, int level)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2, ix3, level);
    if (block) {
        this->blockIdMap.erase(block->getGlobalID());
        return this->levelSet[level].erase(Block3DKey(ix1, ix2, ix3)) > 0;
    } else {
        return false;
    }
}
void Grid3D::deleteBlocks()
{
    std::vector<std::vector<SPtr<Block3D>>> blocksVector(25);
    int minInitLevel = D3Q27System::MINLEVEL;
    int maxInitLevel = D3Q27System::MAXLEVEL;
    for (int level = minInitLevel; level < maxInitLevel; level++) {
        getBlocks(level, blocksVector[level]);
        for (SPtr<Block3D> block : blocksVector[level]) //    blocks of the current level
            deleteBlock(block);
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::replaceBlock(SPtr<Block3D> block)
{
    if (block) {
        this->deleteBlock(block);
        this->addBlock(block);
    }
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::getBlock(int ix1, int ix2, int ix3, int level) const
{
    if (!this->hasLevel(level))
        return SPtr<Block3D>();

    int N1 = (nx1 << level);
    int N2 = (nx2 << level);
    int N3 = (nx3 << level);

    if (!this->isPeriodicX1() && (ix1 > N1 - 1 || ix1 < 0))
        return SPtr<Block3D>();
    else if (this->isPeriodicX1() && (ix1 >= N1 - 1 || ix1 < 0)) {
        ix1 = ((ix1 % N1) + N1) % N1;
    }
    if (!this->isPeriodicX2() && (ix2 > N2 - 1 || ix2 < 0))
        return SPtr<Block3D>();
    else if (this->isPeriodicX2() && (ix2 >= N2 - 1 || ix2 < 0)) {
        ix2 = ((ix2 % N2) + N2) % N2;
    }
    if (!this->isPeriodicX3() && (ix3 > N3 - 1 || ix3 < 0))
        return SPtr<Block3D>();
    else if (this->isPeriodicX3() && (ix3 >= N3 - 1 || ix3 < 0)) {
        ix3 = ((ix3 % N3) + N3) % N3;
    }

    Block3DMap::const_iterator it;
    it = levelSet[level].find(Block3DKey(ix1, ix2, ix3));
    if (it == levelSet[level].end())
        return SPtr<Block3D>();
    else
        return it->second;
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::getBlock(int id) const
{
    BlockIDMap::const_iterator it;
    if ((it = blockIdMap.find(id)) == blockIdMap.end()) {
        return SPtr<Block3D>();
    }

    return it->second;
}
//////////////////////////////////////////////////////////////////////////
Grid3D::BlockIDMap &Grid3D::getBlockIDs() { return blockIdMap; }
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::getSuperBlock(SPtr<Block3D> block)
{
    int ix1   = block->getX1();
    int ix2   = block->getX2();
    int ix3   = block->getX3();
    int level = block->getLevel();
    return getSuperBlock(ix1, ix2, ix3, level);
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::getSuperBlock(int ix1, int ix2, int ix3, int level)
{
    if (!this->hasLevel(level))
        return SPtr<Block3D>();
    if (level < 1)
        throw UbException(UB_EXARGS, "level <1");

    // from Lower Level to higher:     >>     1 in x1,x2,x3
    SPtr<Block3D> block;
    for (int l = level - 1; l >= 0; l--) {
        ix1 = ix1 >> 1;
        ix2 = ix2 >> 1;
        ix3 = ix3 >> 1;

        block = this->getBlock(ix1, ix2, ix3, l);
        if (block)
            return block;
    }
    return SPtr<Block3D>();
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocks(SPtr<Block3D> block, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    int ix1   = block->getX1();
    int ix2   = block->getX2();
    int ix3   = block->getX3();
    int level = block->getLevel();
    getSubBlocks(ix1, ix2, ix3, level, levelDepth, blocks);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocks(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    if (!this->getBlock(ix1, ix2, ix3, level))
        return;
    if (level > 0 && !this->getSuperBlock(ix1, ix2, ix3, level))
        return;
    if (level >= D3Q27System::MAXLEVEL)
        throw UbException(UB_EXARGS, "Level bigger then MAXLEVEL");

    int x1[] = { ix1 << 1, (ix1 << 1) + 1 };
    int x2[] = { ix2 << 1, (ix2 << 1) + 1 };
    int x3[] = { ix3 << 1, (ix3 << 1) + 1 };
    int l    = level + 1;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                SPtr<Block3D> block = this->getBlock(x1[i], x2[j], x3[k], l);
                if (block)
                    blocks.push_back(block);
                else if (l < levelDepth)
                    this->getSubBlocks(x1[i], x2[j], x3[k], l, levelDepth, blocks);
            }
}
//////////////////////////////////////////////////////////////////////////
bool Grid3D::expandBlock(int ix1, int ix2, int ix3, int level)
{
    this->checkLevel(level);

    SPtr<Block3D> block = this->getBlock(ix1, ix2, ix3, level);
    if (!block)
        throw UbException(UB_EXARGS, "block(x1=" + UbSystem::toString(ix1) + ", x2=" + UbSystem::toString(ix2) +
                                         ", x3=" + UbSystem::toString(ix3) + ", l=" + UbSystem::toString(level) +
                                         ") is not exist");

    // da bei periodic der eigentliche block andere indizes hat:
    ix1 = block->getX1();
    ix2 = block->getX2();
    ix3 = block->getX3();

    int l = level + 1;
    if (l > D3Q27System::MAXLEVEL)
        throw UbException(UB_EXARGS, "level > Grid3D::MAXLEVEL");

    int west   = ix1 << 1;
    int east   = west + 1;
    int south  = ix2 << 1;
    int north  = south + 1;
    int bottom = ix3 << 1;
    int top    = bottom + 1;

    auto blockBSW = std::make_shared<Block3D>(west, south, bottom, l);
    auto blockBSE = std::make_shared<Block3D>(east, south, bottom, l);
    auto blockBNW = std::make_shared<Block3D>(west, north, bottom, l);
    auto blockBNE = std::make_shared<Block3D>(east, north, bottom, l);
    auto blockTSW = std::make_shared<Block3D>(west, south, top, l);
    auto blockTSE = std::make_shared<Block3D>(east, south, top, l);
    auto blockTNW = std::make_shared<Block3D>(west, north, top, l);
    auto blockTNE = std::make_shared<Block3D>(east, north, top, l);

    if (!this->deleteBlock(ix1, ix2, ix3, level))
        throw UbException(UB_EXARGS, "could not delete block");

    this->addBlock(blockBSW);
    this->addBlock(blockBSE);
    this->addBlock(blockBNW);
    this->addBlock(blockBNE);
    this->addBlock(blockTSW);
    this->addBlock(blockTSE);
    this->addBlock(blockTNW);
    this->addBlock(blockTNE);

    return true;
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::collapseBlock(int fix1, int fix2, int fix3, int flevel, int levelDepth)
{
    using UbSystem::toString;

    SPtr<Block3D> fblock = this->getBlock(fix1, fix2, fix3, flevel);
    if (flevel < 1)
        throw UbException(UB_EXARGS, "level of block (" + toString(fix1) + "," + toString(fix2) + "," + toString(fix3) +
                                         "," + toString(flevel) + ") is < 1");
    if (!fblock) {
        throw UbException(UB_EXARGS, "specific block(" + toString(fix1) + "," + toString(fix2) + "," + toString(fix3) +
                                         "," + toString(flevel) + ") doesn't exists");
    }
    if (!fblock->isActive())
        throw UbException(UB_EXARGS, "block(" + toString(fix1) + "," + toString(fix2) + "," + toString(fix3) + "," +
                                         toString(flevel) + ") is not active");

    // da bei periodic der eigentliche block andere indizes hat:
    fix1 = fblock->getX1();
    fix2 = fblock->getX2();
    fix3 = fblock->getX3();

    int cix1 = fblock->getX1() >> 1;
    int cix2 = fblock->getX2() >> 1;
    int cix3 = fblock->getX3() >> 1;

    int fx1[2] = { cix1 << 1, (cix1 << 1) + 1 };
    int fx2[2] = { cix2 << 1, (cix2 << 1) + 1 };
    int fx3[2] = { cix3 << 1, (cix3 << 1) + 1 };
    int clevel = flevel - 1;

    vector<SPtr<Block3D>> blocks;
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++) {
                this->getSubBlocks(fx1[k], fx2[i], fx3[l], flevel, levelDepth, blocks);
                while (!blocks.empty()) {
                    // man muss nur eine von den moeglichen acht "collapsen", die anderen werden
                    // dann (rekursiv) collapsed, da die schleife oben alle vier abfragt
                    this->collapseBlock(blocks[0]->getX1(), blocks[0]->getX2(), blocks[0]->getX3(),
                                        blocks[0]->getLevel(), levelDepth);
                    this->getSubBlocks(fx1[k], fx2[i], fx3[l], flevel, levelDepth, blocks);
                }
            }

    vector<SPtr<Block3D>> fineBlocks(8);
    /*BSW*/ fineBlocks[0] = this->getBlock(fx1[0], fx2[0], fx3[0], flevel);
    /*BSE*/ fineBlocks[1] = this->getBlock(fx1[1], fx2[0], fx3[0], flevel);
    /*BNE*/ fineBlocks[2] = this->getBlock(fx1[1], fx2[1], fx3[0], flevel);
    /*BNW*/ fineBlocks[3] = this->getBlock(fx1[0], fx2[1], fx3[0], flevel);
    /*TSW*/ fineBlocks[4] = this->getBlock(fx1[0], fx2[0], fx3[1], flevel);
    /*TSE*/ fineBlocks[5] = this->getBlock(fx1[1], fx2[0], fx3[1], flevel);
    /*TNE*/ fineBlocks[6] = this->getBlock(fx1[1], fx2[1], fx3[1], flevel);
    /*TNW*/ fineBlocks[7] = this->getBlock(fx1[0], fx2[1], fx3[1], flevel);

    auto cblock = std::make_shared<Block3D>(cix1, cix2, cix3, clevel);

    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
                if (!this->deleteBlock(fx1[k], fx2[i], fx3[l], flevel))
                    throw UbException(UB_EXARGS, "could not delete block");

    this->addBlock(cblock);

    return cblock;
}
//////////////////////////////////////////////////////////////////////////
// TODO: make visitor for this
void Grid3D::deleteConnectors()
{
    for (Block3DMap blockMap : levelSet) {
        for (Block3DMap::value_type b : blockMap) {
            SPtr<Block3D> block = b.second;
            block->deleteConnectors();
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::setRank(int rank) { this->rank = rank; }
//////////////////////////////////////////////////////////////////////////
int Grid3D::getRank() const { return rank; }
//////////////////////////////////////////////////////////////////////////
int Grid3D::getBundle() const { return bundle; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setBundle(int bundle) { this->bundle = bundle; }
//////////////////////////////////////////////////////////////////////////
bool Grid3D::isPeriodicX1() const { return this->periodicX1; }
//////////////////////////////////////////////////////////////////////////
bool Grid3D::isPeriodicX2() const { return this->periodicX2; }
//////////////////////////////////////////////////////////////////////////
bool Grid3D::isPeriodicX3() const { return this->periodicX3; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setPeriodicX1(bool value) { this->periodicX1 = value; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setPeriodicX2(bool value) { this->periodicX2 = value; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setPeriodicX3(bool value) { this->periodicX3 = value; }
//////////////////////////////////////////////////////////////////////////
UbTupleInt3 Grid3D::getBlockIndexes(real blockX1Coord, real blockX2Coord, real blockX3Coord) const
{
    if (!trafo) {
        return makeUbTuple((int)blockX1Coord, (int)blockX2Coord, (int)blockX3Coord);
    }

    return makeUbTuple((int)trafo->transformForwardToX1Coordinate(blockX1Coord, blockX2Coord, blockX3Coord),
                       (int)trafo->transformForwardToX2Coordinate(blockX1Coord, blockX2Coord, blockX3Coord),
                       (int)trafo->transformForwardToX3Coordinate(blockX1Coord, blockX2Coord, blockX3Coord));
}
//////////////////////////////////////////////////////////////////////////
UbTupleInt3 Grid3D::getBlockIndexes(real blockX1Coord, real blockX2Coord, real blockX3Coord, int level) const
{
    if (!trafo) {
        return makeUbTuple((int)blockX1Coord, (int)blockX2Coord, (int)blockX3Coord);
    }

    real dx = getDeltaX(level);
    real blockLentghX1, blockLentghX2, blockLentghX3;
    blockLentghX1      = blockNx1 * dx;
    blockLentghX2      = blockNx2 * dx;
    blockLentghX3      = blockNx3 * dx;
    UbTupleDouble3 org = getBlockWorldCoordinates(0, 0, 0, 0);

    SPtr<CoordinateTransformation3D> trafo_temp(new CoordinateTransformation3D(
        val<1>(org), val<2>(org), val<3>(org), blockLentghX1, blockLentghX2, blockLentghX3));

    if (!trafo_temp) {
        return makeUbTuple((int)blockX1Coord, (int)blockX2Coord, (int)blockX3Coord);
    }

    return makeUbTuple((int)trafo_temp->transformForwardToX1Coordinate(blockX1Coord, blockX2Coord, blockX3Coord),
                       (int)trafo_temp->transformForwardToX2Coordinate(blockX1Coord, blockX2Coord, blockX3Coord),
                       (int)trafo_temp->transformForwardToX3Coordinate(blockX1Coord, blockX2Coord, blockX3Coord));
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 Grid3D::getBlockLengths(const SPtr<Block3D> block) const
{
    int level    = block->getLevel();
    real delta = 1.0 / (real)(1 << level);

    if (!trafo)
        makeUbTuple<real, real, real>(delta, delta, delta);

    return makeUbTuple(trafo->getX1CoordinateScaling() * delta, trafo->getX2CoordinateScaling() * delta,
                       trafo->getX3CoordinateScaling() * delta);
}
//////////////////////////////////////////////////////////////////////////
using namespace vf::basics::constant;
UbTupleDouble6 Grid3D::getBlockOversize() const { return makeUbTuple(0., 0., 0., 0., 0., 0.); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setCoordinateTransformator(SPtr<CoordinateTransformation3D> trafo) { this->trafo = trafo; }
//////////////////////////////////////////////////////////////////////////
const SPtr<CoordinateTransformation3D> Grid3D::getCoordinateTransformator() const { return this->trafo; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setDeltaX(real dx) { this->orgDeltaX = dx; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setDeltaX(real worldUnit, real gridUnit) { this->orgDeltaX = worldUnit / gridUnit; }
//////////////////////////////////////////////////////////////////////////
real Grid3D::getDeltaX(int level) const
{
    real delta = this->orgDeltaX / (real)(1 << level);
    return delta;
}
//////////////////////////////////////////////////////////////////////////
real Grid3D::getDeltaX(SPtr<Block3D> block) const { return getDeltaX(block->getLevel()); }
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 Grid3D::getNodeOffset(SPtr<Block3D> block) const
{
    double delta = (double)this->getDeltaX(block);
    return makeUbTuple((double)offset * delta, (double)offset * delta, (double)offset * delta);
}
////////////////////////////////////////////////////////////////////////////
Vector3D Grid3D::getNodeCoordinates(SPtr<Block3D> block, int ix1, int ix2, int ix3) const
{
    UbTupleDouble3 org        = this->getBlockWorldCoordinates(block);
    UbTupleDouble3 nodeOffset = this->getNodeOffset(block);
    real deltaX             = getDeltaX(block);

    real x1 = val<1>(org) - val<1>(nodeOffset) + (real)ix1 * deltaX;
    real x2 = val<2>(org) - val<2>(nodeOffset) + (real)ix2 * deltaX;
    real x3 = val<3>(org) - val<3>(nodeOffset) + (real)ix3 * deltaX;

    return Vector3D(x1, x2, x3);
}
////////////////////////////////////////////////////////////////////////////
UbTupleInt3 Grid3D::getNodeIndexes(SPtr<Block3D> block, real nodeX1Coord, real nodeX2Coord,
                                   real nodeX3Coord) const
{
    UbTupleDouble3 org        = this->getBlockWorldCoordinates(block);
    UbTupleDouble3 nodeOffset = this->getNodeOffset(block);
    real deltaX             = getDeltaX(block);

    int ix1, ix2, ix3;
    real ixx1 = (abs(nodeX1Coord - val<1>(org) + val<1>(nodeOffset)) / deltaX);
    real ixx2 = (abs(nodeX2Coord - val<2>(org) + val<2>(nodeOffset)) / deltaX);
    real ixx3 = (abs(nodeX3Coord - val<3>(org) + val<3>(nodeOffset)) / deltaX);
    if (ixx1 - (int)ixx1 > .9999999999)
        ix1 = (int)ixx1 + 1;
    else
        ix1 = (int)ixx1;
    if (ixx2 - (int)ixx2 > .9999999999)
        ix2 = (int)ixx2 + 1;
    else
        ix2 = (int)ixx2;
    if (ixx3 - (int)ixx3 > .9999999999)
        ix3 = (int)ixx3 + 1;
    else
        ix3 = (int)ixx3;

    return makeUbTuple(ix1, ix2, ix3);
}
//////////////////////////////////////////////////////////////////////////
// returns tuple with origin of block in world-coordinates
UbTupleDouble3 Grid3D::getBlockWorldCoordinates(SPtr<Block3D> block) const
{
    if (!block)
        throw UbException(UB_EXARGS, "block " + block->toString() + "is not exist");

    int blockX1Index = block->getX1();
    int blockX2Index = block->getX2();
    int blockX3Index = block->getX3();
    int level        = block->getLevel();

    return this->getBlockWorldCoordinates(blockX1Index, blockX2Index, blockX3Index, level);
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 Grid3D::getBlockWorldCoordinates(int blockX1Index, int blockX2Index, int blockX3Index, int level) const
{
    real c1oShiftedLevel = 1.0 / (real)(1 << level);
    real x1              = (real)blockX1Index * c1oShiftedLevel;
    real x2              = (real)blockX2Index * c1oShiftedLevel;
    real x3              = (real)blockX3Index * c1oShiftedLevel;

    if (!trafo)
        return { x1, x2, x3 };

    return { trafo->transformBackwardToX1Coordinate(x1, x2, x3), trafo->transformBackwardToX2Coordinate(x1, x2, x3),
             trafo->transformBackwardToX3Coordinate(x1, x2, x3) };
}
//////////////////////////////////////////////////////////////////////////
// double Grid3D::getDeltaT(SPtr<Block3D> block) const
//{
//   int    level = block->getLevel();
//   double delta = 1.0/(double)(1<<level);
//   return delta;
//}
//////////////////////////////////////////////////////////////////////////
void Grid3D::checkLevel(int level)
{
    if (level < 0) {
        throw UbException(UB_EXARGS, "l(" + UbSystem::toString(level) + (string) ")<0");
    }
    if (level > D3Q27System::MAXLEVEL) {
        throw UbException(UB_EXARGS, "l(" + UbSystem::toString(level) + (string) ")>MAXLEVEL");
    }
    if (this->levelSet[level].size() == 0) {
        throw UbException(UB_EXARGS, "levelMap for level(" + UbSystem::toString(level) + (string) ")==NULL");
    }
}
//////////////////////////////////////////////////////////////////////////
bool Grid3D::hasLevel(int level) const
{
    if (level < 0)
        return false;
    if (level > D3Q27System::MAXLEVEL)
        return false;
    if (this->levelSet[level].size() == 0)
        return false;

    return true;
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::setBlockNX(int nx1, int nx2, int nx3)
{
    blockNx1 = nx1;
    blockNx2 = nx2;
    blockNx3 = nx3;
}
//////////////////////////////////////////////////////////////////////////
UbTupleInt3 Grid3D::getBlockNX() const { return makeUbTuple(blockNx1, blockNx2, blockNx3); }
//////////////////////////////////////////////////////////////////////////

SPtr<Block3D> Grid3D::getNeighborBlock(int dir, int ix1, int ix2, int ix3, int level) const
{
    return this->getBlock(ix1 + D3Q27System::DX1[dir], ix2 + D3Q27System::DX2[dir], ix3 + D3Q27System::DX3[dir],
                          level);
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3D> Grid3D::getNeighborBlock(int dir, SPtr<Block3D> block) const
{
    int x1    = block->getX1();
    int x2    = block->getX2();
    int x3    = block->getX3();
    int level = block->getLevel();
    return this->getNeighborBlock(dir, x1, x2, x3, level);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getAllNeighbors(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
   // for (int dir = D3Q27System::STARTDIR; dir <= D3Q27System::ENDDIR; dir++)FSTARTDIR
   for (int dir = D3Q27System::FSTARTDIR; dir <= D3Q27System::FENDDIR; dir++)
   {
        this->getNeighborBlocksForDirection(dir, ix1, ix2, ix3, level, levelDepth, blocks);
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getAllNeighbors(SPtr<Block3D> block, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    int x1 = block->getX1();
    int x2 = block->getX2();
    int x3 = block->getX3();
    getAllNeighbors(x1, x2, x3, level, levelDepth, blocks);
}
//////////////////////////////////////////////////////////////////////////
/**
 * Returns all direct northern neighbor cells of the specified grid cell (may be NULL!).
 * @param ix1 index in x1 direction
 * @param ix2 index in x2 direction
 * @param ix3 index in x3 direction
 * @param level the level
 */
void Grid3D::getNeighborsNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 + 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 + 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksSouth(ix1, ix2 + 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTop(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottom(ix1, ix2, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottom(int ix1, int ix2, int ix3, int level, int levelDepth,
                                std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTop(ix1, ix2, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 - 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 - 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksNorth(ix1, ix2 - 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksWest(ix1 + 1, ix2, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksEast(ix1 - 1, ix2, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//   diagonals
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                   std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 + 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 + 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksSouthWest(ix1 + 1, ix2 + 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                   std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 + 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 + 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksSouthEast(ix1 - 1, ix2 + 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                   std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 - 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 - 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksNorthWest(ix1 + 1, ix2 - 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                   std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 - 1, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 - 1, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksNorthEast(ix1 - 1, ix2 - 1, ix3, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//   diagonals  top
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                 std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomWest(ix1 + 1, ix2, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                 std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomEast(ix1 - 1, ix2, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopNorth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 + 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 + 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomSouth(ix1, ix2 + 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopSouth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 - 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 - 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomNorth(ix1, ix2 - 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//   diagonals  bottom
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                    std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopWest(ix1 + 1, ix2, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                    std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopEast(ix1 - 1, ix2, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomNorth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 + 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 + 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopSouth(ix1, ix2 + 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomSouth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2 - 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2 - 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopNorth(ix1, ix2 - 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                      std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 + 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 + 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomSouthWest(ix1 + 1, ix2 + 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                      std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 + 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 + 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomSouthEast(ix1 - 1, ix2 + 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                      std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 - 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 - 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomNorthWest(ix1 + 1, ix2 - 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsTopSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                      std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 - 1, ix3 + 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 - 1, ix3 + 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksBottomNorthEast(ix1 - 1, ix2 - 1, ix3 + 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                         std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 + 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 + 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopSouthWest(ix1 + 1, ix2 + 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                         std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 + 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 + 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopSouthEast(ix1 - 1, ix2 + 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                         std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 + 1, ix2 - 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 + 1, ix2 - 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopNorthWest(ix1 + 1, ix2 - 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsBottomSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                         std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1 - 1, ix2 - 1, ix3 - 1, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1 - 1, ix2 - 1, ix3 - 1, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    this->getSubBlocksTopNorthEast(ix1 - 1, ix2 - 1, ix3 - 1, level, blocks, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborBlocksForDirection(int dir, int ix1, int ix2, int ix3, int level, int levelDepth,
                                           std::vector<SPtr<Block3D>> &blocks)
{
    using namespace vf::lbm::dir;

    switch (dir) {
        case dP00:
            this->getNeighborsEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM00:
            this->getNeighborsWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0P0:
            this->getNeighborsNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0M0:
            this->getNeighborsSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d00P:
            this->getNeighborsTop(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d00M:
            this->getNeighborsBottom(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPP0:
            this->getNeighborsNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMM0:
            this->getNeighborsSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPM0:
            this->getNeighborsSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMP0:
            this->getNeighborsNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dP0P:
            this->getNeighborsTopEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM0M:
            this->getNeighborsBottomWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dP0M:
            this->getNeighborsBottomEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM0P:
            this->getNeighborsTopWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0PP:
            this->getNeighborsTopNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0MM:
            this->getNeighborsBottomSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0PM:
            this->getNeighborsBottomNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0MP:
            this->getNeighborsTopSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPPP:
            this->getNeighborsTopNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMPP:
            this->getNeighborsTopNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPMP:
            this->getNeighborsTopSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMMP:
            this->getNeighborsTopSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPPM:
            this->getNeighborsBottomNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMPM:
            this->getNeighborsBottomNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPMM:
            this->getNeighborsBottomSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMMM:
            this->getNeighborsBottomSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        default:
            throw UbException(UB_EXARGS, "direction " + UbSystem::toString(dir) + " is not exist");
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborsZero(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks)
{
    SPtr<Block3D> block = this->getBlock(ix1, ix2, ix3, level);
    if (block) {
        blocks.push_back(block);
    }

    if (level > 0) {
        block = this->getSuperBlock(ix1, ix2, ix3, level);
        if (block) {
            blocks.push_back(block);
        }
    }
    // this->getSubBlocksNull(ix1, ix2, ix3, level, blocks, levelDepth);
    this->getSubBlocks(ix1, ix2, ix3, level, levelDepth, blocks);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksZero(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1E = (ix1 << 1) + 1;
    int x1W = (ix1 << 1);
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3B = ix3 << 1;
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1E, x2S, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1E, x2S, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2N, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1E, x2N, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2S, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1E, x2S, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2N, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1E, x2N, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2S, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1W, x2S, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2N, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1W, x2N, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2S, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1W, x2S, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2N, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1W, x2N, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getNeighborBlocksForDirectionWithREST(int dir, int ix1, int ix2, int ix3, int level, int levelDepth,
                                                      std::vector<SPtr<Block3D>> &blocks)
{
    using namespace vf::lbm::dir;

    switch (dir) {
        case dP00:
            this->getNeighborsEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM00:
            this->getNeighborsWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0P0:
            this->getNeighborsNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0M0:
            this->getNeighborsSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d00P:
            this->getNeighborsTop(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d00M:
            this->getNeighborsBottom(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPP0:
            this->getNeighborsNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMM0:
            this->getNeighborsSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPM0:
            this->getNeighborsSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMP0:
            this->getNeighborsNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dP0P:
            this->getNeighborsTopEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM0M:
            this->getNeighborsBottomWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dP0M:
            this->getNeighborsBottomEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dM0P:
            this->getNeighborsTopWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0PP:
            this->getNeighborsTopNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0MM:
            this->getNeighborsBottomSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0PM:
            this->getNeighborsBottomNorth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d0MP:
            this->getNeighborsTopSouth(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPPP:
            this->getNeighborsTopNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMPP:
            this->getNeighborsTopNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPMP:
            this->getNeighborsTopSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMMP:
            this->getNeighborsTopSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPPM:
            this->getNeighborsBottomNorthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMPM:
            this->getNeighborsBottomNorthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dPMM:
            this->getNeighborsBottomSouthEast(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case dMMM:
            this->getNeighborsBottomSouthWest(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        case d000:
            this->getNeighborsZero(ix1, ix2, ix3, level, levelDepth, blocks);
            break;
        default:
            throw UbException(UB_EXARGS, "direction " + UbSystem::toString(dir) + " is not exist");
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksEast(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1  = (ix1 << 1) + 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3B = ix3 << 1;
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1, x2S, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1, x2S, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2N, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1, x2N, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2S, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1, x2S, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2N, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksEast(x1, x2N, x3T, l, blockVector, levelDepth);
}

//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksWest(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1  = ix1 << 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3B = ix3 << 1;
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1, x2S, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksWest(x1, x2S, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2N, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksWest(x1, x2N, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2S, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksWest(x1, x2S, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1, x2N, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksWest(x1, x2N, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksNorth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1W = ix1 << 1;
    int x1E = x1W + 1;
    int x2  = (ix2 << 1) + 1;
    int x3B = ix3 << 1;
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1W, x2, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksNorth(x1W, x2, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksNorth(x1E, x2, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksNorth(x1W, x2, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksNorth(x1E, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksSouth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1W = ix1 << 1;
    int x1E = x1W + 1;
    int x2  = ix2 << 1;
    int x3B = ix3 << 1;
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1W, x2, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksSouth(x1W, x2, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2, x3B, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksSouth(x1E, x2, x3B, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksSouth(x1W, x2, x3T, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2, x3T, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksSouth(x1E, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTop(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector, int levelDepth)
{
    int x1W = ix1 << 1;
    int x1E = x1W + 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3  = (ix3 << 1) + 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1W, x2N, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksTop(x1W, x2N, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2N, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksTop(x1E, x2N, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2S, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksTop(x1W, x2S, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2S, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksTop(x1E, x2S, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottom(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                int levelDepth)
{
    int x1W = ix1 << 1;
    int x1E = x1W + 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3  = ix3 << 1;
    int l   = level + 1;

    SPtr<Block3D> block = this->getBlock(x1W, x2N, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksBottom(x1W, x2N, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2N, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksBottom(x1E, x2N, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1W, x2S, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksBottom(x1W, x2S, x3, l, blockVector, levelDepth);

    block = this->getBlock(x1E, x2S, x3, l);
    if (block != NULL)
        blockVector.push_back(block);
    else if (l < levelDepth)
        this->getSubBlocksBottom(x1E, x2S, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//  diagonals
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksNorthEast(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                   int levelDepth)
{
    int x1  = (ix1 << 1) + 1;
    int x2  = (ix2 << 1) + 1;
    int x3B = (ix3 << 1);
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> blockB = this->getBlock(x1, x2, x3B, l);
    if (blockB)
        blockVector.push_back(blockB);
    else if (l < levelDepth)
        this->getSubBlocksNorthEast(x1, x2, x3B, l, blockVector, levelDepth);

    SPtr<Block3D> blockT = this->getBlock(x1, x2, x3T, l);
    if (blockT)
        blockVector.push_back(blockT);
    else if (l < levelDepth)
        this->getSubBlocksNorthEast(x1, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksNorthWest(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                   int levelDepth)
{
    int x1  = (ix1 << 1);
    int x2  = (ix2 << 1) + 1;
    int x3B = (ix3 << 1);
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> blockB = this->getBlock(x1, x2, x3B, l);
    if (blockB)
        blockVector.push_back(blockB);
    else if (l < levelDepth)
        this->getSubBlocksNorthWest(x1, x2, x3B, l, blockVector, levelDepth);

    SPtr<Block3D> blockT = this->getBlock(x1, x2, x3T, l);
    if (blockT)
        blockVector.push_back(blockT);
    else if (l < levelDepth)
        this->getSubBlocksNorthWest(x1, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksSouthWest(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                   int levelDepth)
{
    int x1  = ix1 << 1;
    int x2  = ix2 << 1;
    int x3B = (ix3 << 1);
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> blockB = this->getBlock(x1, x2, x3B, l);
    if (blockB)
        blockVector.push_back(blockB);
    else if (l < levelDepth)
        this->getSubBlocksSouthWest(x1, x2, x3B, l, blockVector, levelDepth);

    SPtr<Block3D> blockT = this->getBlock(x1, x2, x3T, l);
    if (blockT)
        blockVector.push_back(blockT);
    else if (l < levelDepth)
        this->getSubBlocksSouthWest(x1, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksSouthEast(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                   int levelDepth)
{
    int x1  = (ix1 << 1) + 1;
    int x2  = ix2 << 1;
    int x3B = (ix3 << 1);
    int x3T = x3B + 1;
    int l   = level + 1;

    SPtr<Block3D> blockB = this->getBlock(x1, x2, x3B, l);
    if (blockB)
        blockVector.push_back(blockB);
    else if (l < levelDepth)
        this->getSubBlocksSouthEast(x1, x2, x3B, l, blockVector, levelDepth);

    SPtr<Block3D> blockT = this->getBlock(x1, x2, x3T, l);
    if (blockT)
        blockVector.push_back(blockT);
    else if (l < levelDepth)
        this->getSubBlocksSouthEast(x1, x2, x3T, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//  diagonals
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopEast(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                 int levelDepth)
{
    int x1  = (ix1 << 1) + 1;
    int x2S = (ix2 << 1);
    int x2N = x2S + 1;
    int x3  = (ix3 << 1) + 1;
    int l   = level + 1;

    SPtr<Block3D> blockN = this->getBlock(x1, x2N, x3, l);
    if (blockN)
        blockVector.push_back(blockN);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2N, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockS = this->getBlock(x1, x2S, x3, l);
    if (blockS)
        blockVector.push_back(blockS);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2S, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopWest(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                 int levelDepth)
{
    int x1  = ix1 << 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3  = (ix3 << 1) + 1;
    int l   = level + 1;

    SPtr<Block3D> blockN = this->getBlock(x1, x2N, x3, l);
    if (blockN)
        blockVector.push_back(blockN);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2N, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockS = this->getBlock(x1, x2S, x3, l);
    if (blockS)
        blockVector.push_back(blockS);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2S, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomEast(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                    int levelDepth)
{
    int x1  = (ix1 << 1) + 1;
    int x2S = ix2 << 1;
    int x2N = x2S + 1;
    int x3  = ix3 << 1;
    int l   = level + 1;

    SPtr<Block3D> blockN = this->getBlock(x1, x2N, x3, l);
    if (blockN)
        blockVector.push_back(blockN);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2N, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockS = this->getBlock(x1, x2S, x3, l);
    if (blockS)
        blockVector.push_back(blockS);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2S, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomWest(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                    int levelDepth)
{
    int x1  = (ix1 << 1);
    int x2S = (ix2 << 1);
    int x2N = x2S + 1;
    int x3  = ix3 << 1;
    int l   = level + 1;

    SPtr<Block3D> blockN = this->getBlock(x1, x2N, x3, l);
    if (blockN)
        blockVector.push_back(blockN);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2N, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockS = this->getBlock(x1, x2S, x3, l);
    if (blockS)
        blockVector.push_back(blockS);
    else if (l < levelDepth)
        this->getSubBlocksTopEast(x1, x2S, x3, l, blockVector, levelDepth);
}

//////////////////////////////////////////////////////////////////////////
//  edge-diagonals
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopNorth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth)
{
    int x1E = (ix1 << 1);
    int x1W = x1E + 1;
    int x2  = (ix2 << 1) + 1;
    int x3  = (ix3 << 1) + 1;
    int l   = level + 1;

    SPtr<Block3D> blockE = this->getBlock(x1E, x2, x3, l);
    if (blockE)
        blockVector.push_back(blockE);
    else if (l < levelDepth)
        this->getSubBlocksTopNorth(x1E, x2, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockW = this->getBlock(x1W, x2, x3, l);
    if (blockW)
        blockVector.push_back(blockW);
    else if (l < levelDepth)
        this->getSubBlocksTopNorth(x1W, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopSouth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth)
{
    int x1E = (ix1 << 1);
    int x1W = x1E + 1;
    int x2  = (ix2 << 1);
    int x3  = (ix3 << 1) + 1;
    int l   = level + 1;

    SPtr<Block3D> blockE = this->getBlock(x1E, x2, x3, l);
    if (blockE)
        blockVector.push_back(blockE);
    else if (l < levelDepth)
        this->getSubBlocksTopSouth(x1E, x2, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockW = this->getBlock(x1W, x2, x3, l);
    if (blockW)
        blockVector.push_back(blockW);
    else if (l < levelDepth)
        this->getSubBlocksTopSouth(x1W, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomNorth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth)
{
    int x1E = ix1 << 1;
    int x1W = x1E + 1;
    int x2  = (ix2 << 1) + 1;
    int x3  = ix3 << 1;
    int l   = level + 1;

    SPtr<Block3D> blockE = this->getBlock(x1E, x2, x3, l);
    if (blockE)
        blockVector.push_back(blockE);
    else if (l < levelDepth)
        this->getSubBlocksBottomNorth(x1E, x2, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockW = this->getBlock(x1W, x2, x3, l);
    if (blockW)
        blockVector.push_back(blockW);
    else if (l < levelDepth)
        this->getSubBlocksBottomNorth(x1W, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomSouth(int ix1, int ix2, int ix3, int level, vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth)
{
    int x1E = (ix1 << 1);
    int x1W = x1E + 1;
    int x2  = ix2 << 1;
    int x3  = ix3 << 1;
    int l   = level + 1;

    SPtr<Block3D> blockE = this->getBlock(x1E, x2, x3, l);
    if (blockE)
        blockVector.push_back(blockE);
    else if (l < levelDepth)
        this->getSubBlocksBottomSouth(x1E, x2, x3, l, blockVector, levelDepth);

    SPtr<Block3D> blockW = this->getBlock(x1W, x2, x3, l);
    if (blockW)
        blockVector.push_back(blockW);
    else if (l < levelDepth)
        this->getSubBlocksBottomSouth(x1W, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
//  space-diagonals
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopNorthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                      int levelDepth)
{
    int x1 = (ix1 << 1) + 1;
    int x2 = (ix2 << 1) + 1;
    int x3 = (ix3 << 1) + 1;
    int l  = level + 1;

    SPtr<Block3D> blockTNE = this->getBlock(x1, x2, x3, l);
    if (blockTNE)
        blockVector.push_back(blockTNE);
    else if (l < levelDepth)
        this->getSubBlocksTopNorthEast(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopNorthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                      int levelDepth)
{
    int x1 = ix1 << 1;
    int x2 = (ix2 << 1) + 1;
    int x3 = (ix3 << 1) + 1;
    int l  = level + 1;

    SPtr<Block3D> blockTNW = this->getBlock(x1, x2, x3, l);
    if (blockTNW)
        blockVector.push_back(blockTNW);
    else if (l < levelDepth)
        this->getSubBlocksTopNorthWest(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopSouthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                      int levelDepth)
{
    int x1 = (ix1 << 1) + 1;
    int x2 = ix2 << 1;
    int x3 = (ix3 << 1) + 1;
    int l  = level + 1;

    SPtr<Block3D> blockTNW = this->getBlock(x1, x2, x3, l);
    if (blockTNW)
        blockVector.push_back(blockTNW);
    else if (l < levelDepth)
        this->getSubBlocksTopSouthEast(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksTopSouthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                      int levelDepth)
{
    int x1 = ix1 << 1;
    int x2 = ix2 << 1;
    int x3 = (ix3 << 1) + 1;
    int l  = level + 1;

    SPtr<Block3D> blockTSW = this->getBlock(x1, x2, x3, l);
    if (blockTSW)
        blockVector.push_back(blockTSW);
    else if (l < levelDepth)
        this->getSubBlocksTopSouthWest(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomNorthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                         int levelDepth)
{
    int x1 = (ix1 << 1) + 1;
    int x2 = (ix2 << 1) + 1;
    int x3 = ix3 << 1;
    int l  = level + 1;

    SPtr<Block3D> blockBNE = this->getBlock(x1, x2, x3, l);
    if (blockBNE)
        blockVector.push_back(blockBNE);
    else if (l < levelDepth)
        this->getSubBlocksBottomNorthEast(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomNorthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                         int levelDepth)
{
    int x1 = ix1 << 1;
    int x2 = (ix2 << 1) + 1;
    int x3 = ix3 << 1;
    int l  = level + 1;

    SPtr<Block3D> blockBNW = this->getBlock(x1, x2, x3, l);
    if (blockBNW)
        blockVector.push_back(blockBNW);
    else if (l < levelDepth)
        this->getSubBlocksBottomNorthWest(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomSouthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                         int levelDepth)
{
    int x1 = (ix1 << 1) + 1;
    int x2 = ix2 << 1;
    int x3 = ix3 << 1;
    int l  = level + 1;

    SPtr<Block3D> blockBSE = this->getBlock(x1, x2, x3, l);
    if (blockBSE)
        blockVector.push_back(blockBSE);
    else if (l < levelDepth)
        this->getSubBlocksBottomSouthEast(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getSubBlocksBottomSouthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                         int levelDepth)
{
    int x1 = ix1 << 1;
    int x2 = ix2 << 1;
    int x3 = ix3 << 1;
    int l  = level + 1;

    SPtr<Block3D> blockBSW = this->getBlock(x1, x2, x3, l);
    if (blockBSW)
        blockVector.push_back(blockBSW);
    else if (l < levelDepth)
        this->getSubBlocksBottomSouthWest(x1, x2, x3, l, blockVector, levelDepth);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getBlocks(int level, std::vector<SPtr<Block3D>> &blockVector)
{
    for (Block3DMap::value_type b : levelSet[level]) {
        blockVector.push_back(b.second);
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getBlocks(int level, int rank, std::vector<SPtr<Block3D>> &blockVector)
{
    for (Block3DMap::value_type b : levelSet[level]) {
        SPtr<Block3D> block = b.second;
        int blockRank       = block->getRank();
        if (blockRank == rank) {
            blockVector.push_back(b.second);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getBlocks(int level, int rank, bool active, std::vector<SPtr<Block3D>> &blockVector)
{
    for (Block3DMap::value_type b : levelSet[level]) {
        SPtr<Block3D> block = b.second;
        int blockRank       = block->getRank();

        if (blockRank == rank && active ? block->isActive() : block->isNotActive()) {
            blockVector.push_back(b.second);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
int Grid3D::getFinestInitializedLevel()
{
    for (int i = D3Q27System::MAXLEVEL; i >= 0; i--)
        if (this->levelSet[i].size() > 0)
            return (i);
    return (-1);
}
//////////////////////////////////////////////////////////////////////////
int Grid3D::getCoarsestInitializedLevel()
{
    for (int i = 0; i <= D3Q27System::MAXLEVEL; i++)
        if (this->levelSet[i].size() > 0)
            return (i);
    return (-1);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::setNX1(int nx1) { this->nx1 = nx1; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setNX2(int nx2) { this->nx2 = nx2; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::setNX3(int nx3) { this->nx3 = nx3; }
//////////////////////////////////////////////////////////////////////////
int Grid3D::getNX1() const { return this->nx1; }
//////////////////////////////////////////////////////////////////////////
int Grid3D::getNX2() const { return this->nx2; }
//////////////////////////////////////////////////////////////////////////
int Grid3D::getNX3() const { return this->nx3; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::deleteBlocks(const std::vector<int> &ids)
{
    for (int i : ids) {
        SPtr<Block3D> block = getBlock(i);
        if (block)
            this->deleteBlock(block);
    }
}
//////////////////////////////////////////////////////////////////////////
int Grid3D::getNumberOfBlocks()
{
    int c = 0;
    for (Block3DMap l : levelSet) {
        c += (int)l.size();
    }
    return c;
}
//////////////////////////////////////////////////////////////////////////
int Grid3D::getNumberOfBlocks(int level) { return (int)levelSet[level].size(); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::getBlocksByCuboid(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                               std::vector<SPtr<Block3D>> &blocks)
{
    int coarsestLevel = this->getCoarsestInitializedLevel();
    int finestLevel   = this->getFinestInitializedLevel();

    //////////////////////////////////////////////////////////////////////////
    // MINIMALE BLOCK-INDIZES BESTIMMEN
    //
    // min:
    real dMinX1 = trafo->transformForwardToX1Coordinate(minX1, minX2, minX3) * (1 << finestLevel);
    real dMinX2 = trafo->transformForwardToX2Coordinate(minX1, minX2, minX3) * (1 << finestLevel);
    real dMinX3 = trafo->transformForwardToX3Coordinate(minX1, minX2, minX3) * (1 << finestLevel);

    // Achtung, wenn minX1 genau auf grenze zwischen zwei bloecken -> der "kleinere" muss genommen werden,
    // da beim Transformieren der "groessere" Index rauskommt
    int iMinX1 = (int)dMinX1;
    if (UbMath::zero(dMinX1 - iMinX1))
        iMinX1 -= 1;
    int iMinX2 = (int)dMinX2;
    if (UbMath::zero(dMinX2 - iMinX2))
        iMinX2 -= 1;
    int iMinX3 = (int)dMinX3;
    if (UbMath::zero(dMinX3 - iMinX3))
        iMinX3 -= 1;

    // max (hier kann die Zusatzabfrage vernachlaessigt werden):
    int iMaxX1 = (int)(trafo->transformForwardToX1Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));
    int iMaxX2 = (int)(trafo->transformForwardToX2Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));
    int iMaxX3 = (int)(trafo->transformForwardToX3Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));

    SPtr<Block3D> block;

    // set, um doppelte bloecke zu vermeiden, die u.U. bei periodic auftreten koennen
    std::set<SPtr<Block3D>> blockset;
    for (int level = coarsestLevel; level <= finestLevel; level++) {
        // damit bei negativen werten auch der "kleinere" genommen wird -> floor!
        int minx1 = (int)std::floor((real)iMinX1 / (1 << (finestLevel - level)));
        int minx2 = (int)std::floor((real)iMinX2 / (1 << (finestLevel - level)));
        int minx3 = (int)std::floor((real)iMinX3 / (1 << (finestLevel - level)));

        int maxx1 = iMaxX1 / (1 << (finestLevel - level));
        int maxx2 = iMaxX2 / (1 << (finestLevel - level));
        int maxx3 = iMaxX3 / (1 << (finestLevel - level));

        for (int ix1 = minx1; ix1 <= maxx1; ix1++)
            for (int ix2 = minx2; ix2 <= maxx2; ix2++)
                for (int ix3 = minx3; ix3 <= maxx3; ix3++)
                    if ((block = this->getBlock(ix1, ix2, ix3, level))) {
                        if (block->getRank() == rank) {
                            blockset.insert(block);
                        }
                    }
    }

    blocks.resize(blockset.size());
    std::copy(blockset.begin(), blockset.end(), blocks.begin());
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getBlocksByCuboid(int level, real minX1, real minX2, real minX3, real maxX1, real maxX2,
                               real maxX3, std::vector<SPtr<Block3D>> &blocks)
{
    //////////////////////////////////////////////////////////////////////////
    // MINIMALE BLOCK-INDIZES BESTIMMEN
    //
    // min:
    real dMinX1 = trafo->transformForwardToX1Coordinate(minX1, minX2, minX3) * (1 << level);
    real dMinX2 = trafo->transformForwardToX2Coordinate(minX1, minX2, minX3) * (1 << level);
    real dMinX3 = trafo->transformForwardToX3Coordinate(minX1, minX2, minX3) * (1 << level);

    // Achtung, wenn minX1 genau auf grenze zwischen zwei bloecken -> der "kleinere" muss genommen werden:
    int iMinX1 = (int)dMinX1;
    if (UbMath::zero(dMinX1 - iMinX1))
        iMinX1 -= 1;
    int iMinX2 = (int)dMinX2;
    if (UbMath::zero(dMinX2 - iMinX2))
        iMinX2 -= 1;
    int iMinX3 = (int)dMinX3;
    if (UbMath::zero(dMinX3 - iMinX3))
        iMinX3 -= 1;

    // max:
    int iMaxX1 = (int)(trafo->transformForwardToX1Coordinate(maxX1, maxX2, maxX3) * (1 << level));
    int iMaxX2 = (int)(trafo->transformForwardToX2Coordinate(maxX1, maxX2, maxX3) * (1 << level));
    int iMaxX3 = (int)(trafo->transformForwardToX3Coordinate(maxX1, maxX2, maxX3) * (1 << level));

    // set, um doppelte bloecke zu vermeiden, die u.U. bei periodic auftreten koennen
    std::set<SPtr<Block3D>> blockset;
    SPtr<Block3D> block;

    for (int ix1 = iMinX1; ix1 <= iMaxX1; ix1++)
        for (int ix2 = iMinX2; ix2 <= iMaxX2; ix2++)
            for (int ix3 = iMinX3; ix3 <= iMaxX3; ix3++)
                if ((block = this->getBlock(ix1, ix2, ix3, level))) {
                    if (block->getRank() == rank) {
                        blockset.insert(block);
                    }
                }

    blocks.resize(blockset.size());
    std::copy(blockset.begin(), blockset.end(), blocks.begin());
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::getAllBlocksByCuboid(real minX1, real minX2, real minX3, real maxX1, real maxX2, real maxX3,
                                  std::vector<SPtr<Block3D>> &blocks)
{
    int coarsestLevel = this->getCoarsestInitializedLevel();
    int finestLevel   = this->getFinestInitializedLevel();

    //////////////////////////////////////////////////////////////////////////
    // MINIMALE BLOCK-INDIZES BESTIMMEN
    //
    // min:
    real dMinX1 = trafo->transformForwardToX1Coordinate(minX1, minX2, minX3) * (1 << finestLevel);
    real dMinX2 = trafo->transformForwardToX2Coordinate(minX1, minX2, minX3) * (1 << finestLevel);
    real dMinX3 = trafo->transformForwardToX3Coordinate(minX1, minX2, minX3) * (1 << finestLevel);

    // Achtung, wenn minX1 genau auf grenze zwischen zwei bloecken -> der "kleinere" muss genommen werden,
    // da beim Transformieren der "groessere" Index rauskommt
    int iMinX1 = (int)dMinX1;
    if (UbMath::zero(dMinX1 - iMinX1))
        iMinX1 -= 1;
    int iMinX2 = (int)dMinX2;
    if (UbMath::zero(dMinX2 - iMinX2))
        iMinX2 -= 1;
    int iMinX3 = (int)dMinX3;
    if (UbMath::zero(dMinX3 - iMinX3))
        iMinX3 -= 1;

    // max (hier kann die Zusatzabfrage vernachlaessigt werden):
    int iMaxX1 = (int)(trafo->transformForwardToX1Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));
    int iMaxX2 = (int)(trafo->transformForwardToX2Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));
    int iMaxX3 = (int)(trafo->transformForwardToX3Coordinate(maxX1, maxX2, maxX3) * (1 << finestLevel));

    SPtr<Block3D> block;

    // set, um doppelte bloecke zu vermeiden, die u.U. bei periodic auftreten koennen
    std::set<SPtr<Block3D>> blockset;
    for (int level = coarsestLevel; level <= finestLevel; level++) {
        // damit bei negativen werten auch der "kleinere" genommen wird -> floor!
        int minx1 = (int)std::floor((real)iMinX1 / (1 << (finestLevel - level)));
        int minx2 = (int)std::floor((real)iMinX2 / (1 << (finestLevel - level)));
        int minx3 = (int)std::floor((real)iMinX3 / (1 << (finestLevel - level)));

        int maxx1 = iMaxX1 / (1 << (finestLevel - level));
        int maxx2 = iMaxX2 / (1 << (finestLevel - level));
        int maxx3 = iMaxX3 / (1 << (finestLevel - level));

        for (int ix1 = minx1; ix1 <= maxx1; ix1++)
            for (int ix2 = minx2; ix2 <= maxx2; ix2++)
                for (int ix3 = minx3; ix3 <= maxx3; ix3++)
                    if ((block = this->getBlock(ix1, ix2, ix3, level))) {
                        if (block) {
                            blockset.insert(block);
                        }
                    }
    }

    blocks.resize(blockset.size());
    std::copy(blockset.begin(), blockset.end(), blocks.begin());
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::calcStartCoordinatesAndDelta(SPtr<Block3D> block, real &worldX1, real &worldX2, real &worldX3,
                                          real &deltaX)
{
    int blocklevel = block->getLevel();
    worldX1        = block->getX1() / (float)(1 << blocklevel);
    worldX2        = block->getX2() / (float)(1 << blocklevel);
    worldX3        = block->getX3() / (float)(1 << blocklevel);
    deltaX         = (real)1.0 / (real)(this->blockNx1 * (real)(1 << blocklevel));

    if (this->trafo) {
        real x1tmp = worldX1, x2tmp = worldX2, x3tmp = worldX3;
        worldX1 = this->trafo->transformBackwardToX1Coordinate(x1tmp, x2tmp, x3tmp);
        worldX2 = this->trafo->transformBackwardToX2Coordinate(x1tmp, x2tmp, x3tmp);
        worldX3 = this->trafo->transformBackwardToX3Coordinate(x1tmp, x2tmp, x3tmp);
        deltaX  = this->trafo->getX1CoordinateScaling() / (real)(this->blockNx1 * (real)(1 << blocklevel));
    }
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::calcStartCoordinatesWithOutOverlap(SPtr<Block3D> block, real &worldX1, real &worldX2, real &worldX3)
{
    int blocklevel = block->getLevel();
    worldX1        = block->getX1() / (float)(1 << blocklevel);
    worldX2        = block->getX2() / (float)(1 << blocklevel);
    worldX3        = block->getX3() / (float)(1 << blocklevel);

    if (this->trafo) {
        real x1tmp = worldX1, x2tmp = worldX2, x3tmp = worldX3;
        worldX1 = this->trafo->transformBackwardToX1Coordinate(x1tmp, x2tmp, x3tmp);
        worldX2 = this->trafo->transformBackwardToX2Coordinate(x1tmp, x2tmp, x3tmp);
        worldX3 = this->trafo->transformBackwardToX3Coordinate(x1tmp, x2tmp, x3tmp);
    }
}
//////////////////////////////////////////////////////////////////////////
int Grid3D::getGhostLayerWidth() const
{
    return static_cast<int>(offset + 0.5);
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::setGhostLayerWidth(int ghostLayerWidth)
{
    this->offset = static_cast<real>(ghostLayerWidth) - c1o2;
}
//////////////////////////////////////////////////////////////////////////
void Grid3D::setTimeStep(real step) { timeStep = step; }
//////////////////////////////////////////////////////////////////////////
real Grid3D::getTimeStep() const { return timeStep; }
//////////////////////////////////////////////////////////////////////////
void Grid3D::fillExtentWithBlocks(UbTupleInt3 minInd, UbTupleInt3 maxInd)
{
    for (int x3 = val<3>(minInd); x3 < val<3>(maxInd); x3++) {
        for (int x2 = val<2>(minInd); x2 < val<2>(maxInd); x2++) {
            for (int x1 = val<1>(minInd); x1 < val<1>(maxInd); x1++) {
                SPtr<Block3D> block(new Block3D(x1, x2, x3, 0));
                this->addBlock(block);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
// void Grid3D::notifyObservers( double step )
//{
//   for(ObserverPtr o, observers)
//   {
//      o->update(step);
//   }
//
//   //std::list<ObserverWeakPtr>::iterator iter = observers.begin();
//
//   //GridObserversSet::iterator iter = observers.begin();
//   //while(iter != observers.end())
//   //{
//   //   if ((*iter).expired())
//   //   {
//   //      iter = observers.erase(iter);
//   //   }
//   //   else
//   //   {
//   //      ObserverPtr observer = (*iter).lock(); // create a shared_ptr from the weak_ptr
//   //      observer->update(step);
//   //      ++iter;
//   //   }
//   //}
//
//}
//////////////////////////////////////////////////////////////////////////
// void Grid3D::addObserver( ObserverPtr observer )
//{
//   observers.insert(observer);
//   //observers.push_back(observer);
//}
////////////////////////////////////////////////////////////////////////////
// void Grid3D::removeObserver( ObserverPtr observer )
//{
//   observers.erase(observer);
//   //observers.remove(observer);
//}
//////////////////////////////////////////////////////////////////////////
void Grid3D::deleteBlockIDs() { this->blockIdMap.clear(); }
//////////////////////////////////////////////////////////////////////////
void Grid3D::renumberBlockIDs()
{
    deleteBlockIDs();

    int startLevel = getCoarsestInitializedLevel();
    int stopLevel  = getFinestInitializedLevel();
    int counter    = 0;

    for (int l = startLevel; l <= stopLevel; l++) {
        std::vector<SPtr<Block3D>> blockVector;
        getBlocks(l, blockVector);
        for (SPtr<Block3D> block : blockVector) {
            block->setGlobalID(counter);
            blockIdMap.insert(std::make_pair(counter, block));
            // Block3D::setMaxGlobalID(counter);
            counter++;
        }
    }
}


//////////////////////////////////////////////////////////////////////////
void Grid3D::updateDistributedBlocks(std::shared_ptr<vf::parallel::Communicator> comm)
{

    std::vector<int> blocks;

    if (comm->isRoot()) {
        int startLevel = getCoarsestInitializedLevel();
        int stopLevel  = getFinestInitializedLevel();

        for (int l = startLevel; l <= stopLevel; l++) {
            std::vector<SPtr<Block3D>> blockVector;
            getBlocks(l, blockVector);
            for (SPtr<Block3D> block : blockVector) {
                blocks.push_back(block->getX1());
                blocks.push_back(block->getX2());
                blocks.push_back(block->getX3());
                blocks.push_back(l);
                blocks.push_back(block->getGlobalID());
            }
        }
    }

    comm->broadcast(blocks);

    if (!comm->isRoot()) {
        int startLevel = getCoarsestInitializedLevel();
        int stopLevel  = getFinestInitializedLevel();

        blockIdMap.clear();

        for (int l = startLevel; l <= stopLevel; l++) {
            levelSet[l].clear();
        }
        this->levelSet.clear();
        levelSet.resize(D3Q27System::MAXLEVEL + 1);

        int rsize = (int)blocks.size();
        for (int i = 0; i < rsize; i += 5) {
            SPtr<Block3D> block(new Block3D(blocks[i], blocks[i + 1], blocks[i + 2], blocks[i + 3]));
            block->setGlobalID(blocks[i + 4]);
            this->addBlock(block);
        }
    }
}

//////////////////////////////////////////////////////////////////////////
