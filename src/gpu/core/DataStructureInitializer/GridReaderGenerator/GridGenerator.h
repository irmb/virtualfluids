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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef GridReaderGenerator_H
#define GridReaderGenerator_H

#include "../GridProvider.h"

#include <vector>
#include <string>

#include "Calculation/Calculation.h"


class Parameter;
class GridBuilder;
class IndexRearrangementForStreams;
class InterpolationCellGrouper;
namespace vf::parallel
{
class Communicator;
}

//! \class GridGenerator derived class of GridProvider
//! \brief mapping the grid of grid generator to data structure for simulation
class GridGenerator
    : public GridProvider
{
private:
    //! \brief string vector with channel direction
    std::vector<std::string> channelDirections;
    //! \brief string vector with channel direction (boundary conditions)
    std::vector<std::string> channelBoundaryConditions;

    std::shared_ptr<GridBuilder> builder;
    std::unique_ptr<const IndexRearrangementForStreams> indexRearrangement;
    std::unique_ptr<const InterpolationCellGrouper> interpolationGrouper;
    const uint mpiProcessID;

public:
    GridGenerator(std::shared_ptr<GridBuilder> builder, std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager, vf::parallel::Communicator& communicator);
    ~GridGenerator() override;
    //! \brief overwrites the default IndexRearrangementForStreams
    void setIndexRearrangementForStreams(std::unique_ptr<IndexRearrangementForStreams>&& indexRearrangement);

    //! \brief allocates and initialized the data structures for Coordinates and node types
    void allocArrays_CoordNeighborGeo() override;
    //! \brief allocates and initialized the values at the boundary conditions
    void allocArrays_BoundaryValues() override;
    //! \brief allocates and initialized the sub-grid distances at the boundary conditions
    void allocArrays_BoundaryQs() override;
    void allocArrays_OffsetScale() override;
    void allocArrays_taggedFluidNodes() override;

    void tagFluidNodeIndices(const std::vector<uint>& taggedFluidNodeIndices, CollisionTemplate tag, uint level) override;
    void sortFluidNodeTags() override;

    virtual void setDimensions() override;
    virtual void setBoundingBox() override;

    virtual void initPeriodicNeigh(std::vector<std::vector<std::vector<unsigned int> > > periodV, std::vector<std::vector<unsigned int> > periodIndex, std::string way) override;
    
private:
    void setPressureValues(int channelSide) const;
    void setPressRhoBC(int sizePerLevel, int level, int channelSide) const;

    void setVelocityValues(int channelSide) const;
    void setVelocity(int level, int sizePerLevel, int channelSide) const;

    void setOutflowValues(int channelSide) const;
    void setOutflow(int level, int sizePerLevel, int channelSide) const;

    void setPressQs(int channelSide) const;
    void setVelocityQs(int channelSide) const;
    void setOutflowQs(int channelSide) const;
    void setNoSlipQs(int channelSide) const;
    void setGeoQs() const;
    void modifyQElement(int channelSide, unsigned int level) const;

    void initalQStruct(QforBoundaryConditions& Q,int channelSide, unsigned int level) const;
    void printQSize(std::string bc,int channelSide, unsigned int level) const;
    void setSizeNoSlip(int channelSide, unsigned int level) const;
    void setSizeGeoQs(unsigned int level) const;
    void setQ27Size(QforBoundaryConditions &Q, real* QQ, unsigned int sizeQ) const;
    bool hasQs(int channelSide, unsigned int level) const;

    void initalValuesDomainDecompostion();
public:
    void initalGridInformations() override;


private:
    //! \brief verifies if there are invalid nodes, stopper nodes or wrong neighbors
    std::string verifyNeighborIndices(int level) const;
    //! \brief verifies single neighbor index
    //! \param index type integer
    //! \param invalidNodes reference to invalid nodes
    //! \param stopperNodes reference to stopper nodes
    //! \param wrongNeighbors reference to wrong neighbors
    std::string verifyNeighborIndex(int level, int index , int &invalidNodes, int &stopperNodes, int &wrongNeighbors) const;
    //! \brief check the neighbors
    //! \param x,y,z lattice node position
    //! \param numberOfWrongNeighbors reference to the number of wrong neighbors
    //! \param neighborIndex index of neighbor node
    //! \param neighborX,neighborY,neighborZ neighbor lattice node position
    //! \param direction type string
    std::string checkNeighbor(int level, real x, real y, real z, int index, int& numberOfWrongNeihgbors, int neighborIndex, real neighborX, real neighborY, real neighborZ, std::string direction) const;
    //! \brief create the pointers in the struct for the BoundaryConditions from the boundary condition array
    //! \param boundaryConditionStruct is a struct containing information about the boundary condition
    //! \param subgridDistances is a pointer to an array containing the subgrid distances
    //! \param numberOfBCnodes is the number of lattice nodes in the boundary condition
    static void getPointersToBoundaryConditions(QforBoundaryConditions& boundaryConditionStruct, real* subgridDistances, const unsigned int numberOfBCnodes);

private:
    friend class GridGeneratorTests_initalValuesDomainDecompostion;
};

#endif
