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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef MetisPartitioningGridVisitor_h
#define MetisPartitioningGridVisitor_h

#if defined VF_METIS && defined VF_MPI

#include <PointerDefinitions.h>
#include <vector>

#include "Grid3DVisitor.h"
#include "MetisPartitioner.h"

namespace vf::parallel {class Communicator;}

////////////////////////////////////////////////////////////////////////
//! \brief The class implements domain decomposition with METIS library
//! \author Kostyantyn Kucher
//////////////////////////////////////////////////////////////////////////

class Grid3D;
class Block3D;

class MetisPartitioningGridVisitor : public Grid3DVisitor
{
public:
    //! This describe different types of decomposition
    enum GraphType { LevelIntersected, LevelBased };

public:
    //! Constructor
    //! \param comm - communicator
    //! \param graphType - type of decomposition
    //! \param numOfDirs - maximum number of neighbors for each process
    //! \param threads - on/off decomposition for threads
    //! \param numberOfThreads - number of threads
    MetisPartitioningGridVisitor(std::shared_ptr<vf::parallel::Communicator> comm, GraphType graphType, int numOfDirs,
                                 MetisPartitioner::PartType partType = MetisPartitioner::KWAY, bool threads = false,
                                 int numberOfThreads = 0);
    ~MetisPartitioningGridVisitor() override;
    void visit(SPtr<Grid3D> grid) override;
    void setNumberOfProcesses(int np);

protected:
    enum PartLevel { BUNDLE, PROCESS, THREAD };
    void collectData(SPtr<Grid3D> grid, int nofSegments, PartLevel level);
    void buildMetisGraphLevelIntersected(SPtr<Grid3D> grid, int nofSegments, PartLevel level);
    void buildMetisGraphLevelBased(SPtr<Grid3D> grid, int nofSegments, PartLevel level);
    bool getPartitionCondition(SPtr<Block3D> block, PartLevel level);
    void distributePartitionData(SPtr<Grid3D> grid, PartLevel level);
    void clear();
    int getEdgeWeight(int dir);
    int nofSegments;
    int numOfDirs;
    std::vector<int> blockID;
    std::vector<idx_t> parts;
    std::shared_ptr<vf::parallel::Communicator> comm;
    int bundleRoot;
    int processRoot;
    int bundleID;
    int processID;
    int numberOfBundles;
    int numberOfThreads;
    bool threads;
    GraphType graphType;
    MetisPartitioner::PartType partType;
    int numberOfProcesses;
};

#endif // VF_MPI
#endif

//! \}
