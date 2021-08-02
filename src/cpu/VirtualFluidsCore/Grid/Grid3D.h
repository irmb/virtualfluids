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
//! \file Grid.h
//! \ingroup Grid
//! \author Konstantin Kutscher, Soeren Freudiger, Sebastian Geller
//=======================================================================================

#ifndef GRID3D_H
#define GRID3D_H

#include <PointerDefinitions.h>
#include <map>
#include <vector>

#include <basics/utilities/UbKeys.h>
#include <basics/utilities/UbTuple.h>
#include <basics/utilities/Vector3D.h>

class CoordinateTransformation3D;

#include <Block3DVisitor.h>
#include <Grid3DVisitor.h>

class Communicator;
class Block3D;
class Interactor3D;

#define OFFSET 1.5

//! A class implements block grid
//////////////////////////////////////////////////////////////////////////
class Grid3D : public enableSharedFromThis<Grid3D>
{
public:
    using Block3DKey      = UbKeys::Key3<int>;
    using Block3DMap      = std::map<Block3DKey, SPtr<Block3D>>;
    using BlockIDMap      = std::map<int, SPtr<Block3D>>;
    using LevelSet        = std::vector<Block3DMap>;
    using Interactor3DSet = std::vector<SPtr<Interactor3D>>;

public:
    Grid3D();
    Grid3D(SPtr<Communicator> comm);
    Grid3D(SPtr<Communicator> comm, int blockNx1, int blockNx2, int blockNx3, int gridNx1, int gridNx2, int gridNx3);
    virtual ~Grid3D() = default;
    //////////////////////////////////////////////////////////////////////////
    // blocks control
    void addBlock(SPtr<Block3D> block);
    bool deleteBlock(SPtr<Block3D> block);
    bool deleteBlock(int ix1, int ix2, int ix3, int level);
    void deleteBlocks();
    void deleteBlocks(const std::vector<int> &ids);
    void replaceBlock(SPtr<Block3D> block);
    SPtr<Block3D> getBlock(int ix1, int ix2, int ix3, int level) const;
    SPtr<Block3D> getBlock(int id) const;
    void getBlocksByCuboid(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3,
                           std::vector<SPtr<Block3D>> &blocks);
    void getBlocksByCuboid(int level, double minX1, double minX2, double minX3, double maxX1, double maxX2,
                           double maxX3, std::vector<SPtr<Block3D>> &blocks);
    void getAllBlocksByCuboid(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3,
                              std::vector<SPtr<Block3D>> &blocks);
    //! get blocks for level
    void getBlocks(int level, std::vector<SPtr<Block3D>> &blockVector);
    //! get blocks for level with current rank
    void getBlocks(int level, int rank, std::vector<SPtr<Block3D>> &blockVector);
    //! get only active or not active blocks
    void getBlocks(int level, int rank, bool active, std::vector<SPtr<Block3D>> &blockVector);
    int getNumberOfBlocks();
    int getNumberOfBlocks(int level);
    BlockIDMap &getBlockIDs();
    void deleteBlockIDs();
    void renumberBlockIDs();
    void updateDistributedBlocks(SPtr<Communicator> comm);
    SPtr<Block3D> getSuperBlock(SPtr<Block3D> block);
    SPtr<Block3D> getSuperBlock(int ix1, int ix2, int ix3, int level);
    void getSubBlocks(SPtr<Block3D> block, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getSubBlocks(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blockVector);
    SPtr<Block3D> getNeighborBlock(int dir, int ix1, int ix2, int ix3, int level) const;
    SPtr<Block3D> getNeighborBlock(int dir, SPtr<Block3D> block) const;
    bool expandBlock(int ix1, int ix2, int ix3, int level);
    SPtr<Block3D> collapseBlock(int fix1, int fix2, int fix3, int flevel, int levelDepth);
    void getAllNeighbors(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getAllNeighbors(SPtr<Block3D> block, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborBlocksForDirection(int dir, int ix1, int ix2, int ix3, int level, int levelDepth,
                                       std::vector<SPtr<Block3D>> &blocks);
    void getNeighborBlocksForDirectionWithDirZero(int dir, int ix1, int ix2, int ix3, int level, int levelDepth,
                                                  std::vector<SPtr<Block3D>> &blocks);

    void getNeighborsZero(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTop(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottom(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);

    void getNeighborsNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                               std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                               std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                               std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                               std::vector<SPtr<Block3D>> &blocks);

    void getNeighborsTopNorth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopSouth(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopEast(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopWest(int ix1, int ix2, int ix3, int level, int levelDepth, std::vector<SPtr<Block3D>> &blocks);

    void getNeighborsBottomNorth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                 std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomSouth(int ix1, int ix2, int ix3, int level, int levelDepth,
                                 std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                std::vector<SPtr<Block3D>> &blocks);

    void getNeighborsTopNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsTopSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                  std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomNorthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomNorthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomSouthEast(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks);
    void getNeighborsBottomSouthWest(int ix1, int ix2, int ix3, int level, int levelDepth,
                                     std::vector<SPtr<Block3D>> &blocks);
    //////////////////////////////////////////////////////////////////////////
    // level control
    int getFinestInitializedLevel();
    int getCoarsestInitializedLevel();
    //////////////////////////////////////////////////////////////////////////
    void deleteConnectors();
    //////////////////////////////////////////////////////////////////////////
    // interactors control
    void addInteractor(SPtr<Interactor3D> interactor);
    void addAndInitInteractor(SPtr<Interactor3D> interactor, double timestep = 0);
    Interactor3DSet getInteractors();
    //////////////////////////////////////////////////////////////////////////
    // visitors
    void accept(Block3DVisitor &blockVisitor);
    void accept(Grid3DVisitor &gridVisitor);
    void accept(SPtr<Grid3DVisitor> gridVisitor);
    //////////////////////////////////////////////////////////////////////////
    // bundle and rank for distributed memory
    void setBundle(int bundle);
    int getBundle() const;
    int getRank() const;
    void setRank(int rank);
    //////////////////////////////////////////////////////////////////////////
    // periodic boundary
    bool isPeriodicX1() const;
    bool isPeriodicX2() const;
    bool isPeriodicX3() const;
    void setPeriodicX1(bool value);
    void setPeriodicX2(bool value);
    void setPeriodicX3(bool value);
    //////////////////////////////////////////////////////////////////////////
    // Topology
    UbTupleInt3 getBlockIndexes(double blockX1Coord, double blockX2Coord, double blockX3Coord) const;
    UbTupleInt3 getBlockIndexes(double blockX1Coord, double blockX2Coord, double blockX3Coord, int level) const;
    UbTupleDouble3 getBlockLengths(SPtr<Block3D> block) const;
    UbTupleDouble6 getBlockOversize() const;
    void setCoordinateTransformator(SPtr<CoordinateTransformation3D> trafo);
    const SPtr<CoordinateTransformation3D> getCoordinateTransformator() const;
    void setDeltaX(double dx);
    void setDeltaX(double worldUnit, double gridUnit);
    double getDeltaX(int level) const;
    double getDeltaX(SPtr<Block3D> block) const;
    UbTupleDouble3 getNodeOffset(SPtr<Block3D> block) const;
    Vector3D getNodeCoordinates(SPtr<Block3D> block, int ix1, int ix2, int ix3) const;
    UbTupleInt3 getNodeIndexes(SPtr<Block3D> block, double nodeX1Coord, double nodeX2Coord, double nodeX3Coord) const;
    void setBlockNX(int nx1, int nx2, int nx3);
    UbTupleInt3 getBlockNX() const;
    UbTupleDouble3 getBlockWorldCoordinates(SPtr<Block3D> block) const;
    UbTupleDouble3 getBlockWorldCoordinates(int blockX1Index, int blockX2Index, int blockX3Index, int level) const;
    void setNX1(int nx1);
    void setNX2(int nx2);
    void setNX3(int nx3);
    int getNX1() const;
    int getNX2() const;
    int getNX3() const;
    void calcStartCoordinatesAndDelta(SPtr<Block3D> block, double &worldX1, double &worldX2, double &worldX3,
                                      double &deltaX);
    void calcStartCoordinatesWithOutOverlap(SPtr<Block3D> block, double &worldX1, double &worldX2, double &worldX3);
    //////////////////////////////////////////////////////////////////////////
    // LBM
    // double getDeltaT(SPtr<Block3D>) const;
    //////////////////////////////////////////////////////////////////////////
    void setTimeStep(double step);
    double getTimeStep() const;

protected:
    void checkLevel(int level);
    bool hasLevel(int level) const;

    void fillExtentWithBlocks(UbTupleInt3 minInd, UbTupleInt3 maxInd);

    void getSubBlocksZero(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                          int levelDepth);

    void getSubBlocksEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                          int levelDepth);
    void getSubBlocksWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                          int levelDepth);
    void getSubBlocksNorth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                           int levelDepth);
    void getSubBlocksSouth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                           int levelDepth);
    void getSubBlocksTop(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector, int levelDepth);
    void getSubBlocksBottom(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                            int levelDepth);

    void getSubBlocksSouthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                               int levelDepth);
    void getSubBlocksSouthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                               int levelDepth);
    void getSubBlocksNorthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                               int levelDepth);
    void getSubBlocksNorthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                               int levelDepth);

    void getSubBlocksTopEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                             int levelDepth);
    void getSubBlocksTopWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                             int levelDepth);
    void getSubBlocksTopNorth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                              int levelDepth);
    void getSubBlocksTopSouth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                              int levelDepth);

    void getSubBlocksBottomEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                int levelDepth);
    void getSubBlocksBottomWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                int levelDepth);
    void getSubBlocksBottomNorth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                 int levelDepth);
    void getSubBlocksBottomSouth(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                 int levelDepth);

    void getSubBlocksTopNorthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth);
    void getSubBlocksTopNorthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth);
    void getSubBlocksTopSouthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth);
    void getSubBlocksTopSouthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                  int levelDepth);
    void getSubBlocksBottomNorthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth);
    void getSubBlocksBottomNorthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth);
    void getSubBlocksBottomSouthEast(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth);
    void getSubBlocksBottomSouthWest(int ix1, int ix2, int ix3, int level, std::vector<SPtr<Block3D>> &blockVector,
                                     int levelDepth);

private:
    LevelSet levelSet;
    BlockIDMap blockIdMap;
    Interactor3DSet interactors;

    int rank{ 0 };
    int bundle{ 0 };

    bool periodicX1{ false };
    bool periodicX2{ false };
    bool periodicX3{ false };

    int blockNx1{ 0 };
    int blockNx2{ 0 };
    int blockNx3{ 0 };

    int nx1{ 0 };
    int nx2{ 0 };
    int nx3{ 0 };

    SPtr<CoordinateTransformation3D> trafo;
    double orgDeltaX{ 1.0 };

    double timeStep{ 0.0 };
};

#endif
