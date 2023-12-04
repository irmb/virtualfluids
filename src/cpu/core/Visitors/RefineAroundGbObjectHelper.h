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
//! \file RefineAroundGbObjectHelper.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RefineAroundGbObjectHelper_H
#define RefineAroundGbObjectHelper_H

#include <PointerDefinitions.h>
#include "lbm/constants/D3Q27.h"

class Grid3D;
namespace vf::parallel {class Communicator;}
class D3Q27TriFaceMeshInteractor;

//! \brief Refine blocks on base of bounding boxes.
//! \details You need to use <i>addGbObject()</i> to add corresponding bounding boxes. Then call <i>refine()</i>.
//! \author K. Kucher
class RefineAroundGbObjectHelper
{
public:
    //! Constructor
    //! \param grid a smart pointer to the grid object
    //! \param maxRefineLevel an integer for maximal refinement level
    //! \param objectIter a D3Q27TriFaceMeshInteractor object - represent geometry which should be refinement
    //! \param startDistance start distance from geometry for refinement
    //! \param stopDistance stop distance from geometry for refinement
    RefineAroundGbObjectHelper(SPtr<Grid3D> grid, int maxRefineLevel, SPtr<D3Q27TriFaceMeshInteractor> objectIter,
                               real startDistance, real stopDistance, std::shared_ptr<vf::parallel::Communicator> comm);
    virtual ~RefineAroundGbObjectHelper();
    //! start refinement
    void refine();

private:
    SPtr<Grid3D> grid;
    SPtr<D3Q27TriFaceMeshInteractor> objectIter;
    int refineLevel;
    real startDistance, stopDistance;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif
