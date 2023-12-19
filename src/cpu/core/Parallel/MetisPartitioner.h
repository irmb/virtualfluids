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
//! \addtogroup cpu_Parallel Parallel
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================


#ifndef METISPARTITIONER_H
#define METISPARTITIONER_H

#if defined VF_METIS

#include "basics/utilities/UbSystem.h"
#include "metis.h"
#include <PointerDefinitions.h>
#include <string>
#include <vector>


//! \brief Class use METIS library for graph-based partitioning.
//! \author Konstantin Kutscher
//! \date 22.09.2011


class MetisPartitioner
{
public:
    enum PartType { RECURSIVE, KWAY };

public:
    MetisPartitioner();
    idx_t *getMetisOptions();
    void setMetisOptions(int option, idx_t value);
    int partition(int nofParts, PartType ptype);

public:
    std::vector<idx_t> xadj; // adjncy offset of nodes
    //(size = n+1, n=nofNodes)
    std::vector<idx_t> adjncy; // array that stores the adjacency lists of nodes
    //(size = m*2, m= nofEdged, factor 2 because edge A->B AND B->A has to be stored)
    std::vector<idx_t> vwgt;   // vertex weights (size=n*ncon, ncon=nofNodeWeightsPerNode)
    std::vector<idx_t> adjwgt; // array that stores the weights of the adjacency lists
    // (size=2*m)
    idx_t *vsize; // array that stores the computation weights per node
    // (size=n)

    real_t *tpwgts; // holds the wished fratcion of segment i, e.g. tpwgts={0.2, 0.2, 0.6}
    // -> patisions 0 and one will get 20% of the weight each and part 3 60%!
    // (size=nofPartitions)  sum of tpwgts must be 1.0

    real_t *
        ubvec; // This is an array of size ncon that specifies the allowed load imbalance tolerance for each constraint.
    // For the ith partition and jth constraint the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction
    // of the jths constraint total weight. The load imbalances must be greater than 1.0.
    // A NULL value can be passed indicating that the load imbalance tolerance for each constraint should
    // be 1.001 (for ncon=1) or 1.01 (for ncon1).

    std::vector<idx_t>
        part; // This is a vector of size n that upon successful completion stores the partition vector of the graph.
              // The numbering of this vector starts from 0
private:
    idx_t options[METIS_NOPTIONS];
};

#endif

#endif

//! \}
