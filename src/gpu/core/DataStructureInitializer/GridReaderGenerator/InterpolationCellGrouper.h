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
//! \author Anna Wellmann
//! \details See [master thesis of Anna Wellmann]
//=======================================================================================

#ifndef InterpolationCellGrouper_H
#define InterpolationCellGrouper_H

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>
#include <memory>
#include <vector>

struct LBMSimulationParameter;
class GridBuilder;

using LBMSimulationParameters = std::vector<std::shared_ptr<LBMSimulationParameter>>;

class InterpolationCellGrouper {
public:
    //! \brief Construct InterpolationCellGrouper object
    InterpolationCellGrouper(const LBMSimulationParameters &parHs, const LBMSimulationParameters &parDs,
                             SPtr<GridBuilder> builder);

    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief Split the interpolation cells from coarse to fine into border an bulk
    //! \details For communication hiding, the interpolation cells from the coarse to the fine grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! see [master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")]
    void splitCoarseToFineIntoBorderAndBulk(uint level) const;

    //! \brief Split the interpolation cells from fine to coarse into border an bulk
    //! \details For communication hiding, the interpolation cells from the fine to the coarse grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! See [master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")]
    void splitFineToCoarseIntoBorderAndBulk(uint level) const;

protected:
    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief This function reorders the arrays of CFC/CFF indices and sets the pointers and sizes of the new
    //! subarrays: \details The coarse cells for interpolation from coarse to fine (coarseToFineCoarse) are divided into two
    //! subgroups: border and bulk. The fine cells (coarseToFineFine) are reordered accordingly. The offset cells (xOffCF,
    //! yOffCF, zOffCF) must be reordered in the same way.
    void reorderCoarseToFineIntoBorderAndBulk(uint level) const;

    //! \brief This function reorders the arrays of FCC/FCF indices and return pointers and sizes of the new subarrays:
    //! \details The coarse cells for interpolation from fine to coarse (fineToCoarseCoarse) are divided into two subgroups:
    //! border and bulk. The fine cells (fineToCoarseFine) are reordered accordingly. The offset cells (xOffFC,
    //! yOffFC, zOffFC) must be reordered in the same way.
    void reorderFineToCoarseIntoBorderAndBulk(uint level) const;

private:
    SPtr<GridBuilder> builder;
    const LBMSimulationParameters &parHs;
    const LBMSimulationParameters &parDs;
};

#endif
