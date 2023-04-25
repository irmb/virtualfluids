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
//! \file CoProcessor.h
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef CoProcessor_H
#define CoProcessor_H

#include <PointerDefinitions.h>
#include "lbm/constants/D3Q27.h"

class Grid3D;
class UbScheduler;

//! \class CoProcessor
//! \brief An abstract class implements observer design pettern
class CoProcessor
{
public:
    //! Class default constructor
    CoProcessor();
    //! \brief Construct CoProcessor object for grid object and scheduler object.
    //! \pre The Grid3D and UbScheduler objects must exist.
    //! \param grid is observable Grid3D object
    //! \param s is UbScheduler object for scheduling of observer
    //! \details
    //! Class CoProcessor implements the observer design pettern. CoProcessor object is observer.  Grid3D object is
    //! observable.
    CoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
    //! Class destructor
    virtual ~CoProcessor();
    //! \brief Updates observer
    //! \param step is the actual time step
    virtual void process(real step) = 0;

protected:
    SPtr<Grid3D> grid;
    SPtr<UbScheduler> scheduler;
};
#endif
