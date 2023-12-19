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

#ifndef CREATETRANSMITTERSHELPER_H
#define CREATETRANSMITTERSHELPER_H

#include "Block3D.h"
#include <parallel/Communicator.h>

#include "LBMSystem.h"

#include <basics/container/CbVector.h>

#include <parallel/transmitter/TbTransmitter.h>
#include <parallel/transmitter/TbTransmitterMpiPool.h>

//! \brief The class helps to create Transmitters.
//! \details It is created two types of Transmitters: MPI and BOND
class CreateTransmittersHelper
{
public:
    //! Switch between same level and interpolation Connectors. NONE - for same level Connectors; SW, NW, NE, SE -
    //! source/target fine blocks in grid interface
    enum IBlock { NONE, SW, NW, NE, SE };
    //! Switch between MPI and BOND Transmitters
    enum TransmitterType { MPI, BOND, MPI2BOND };

public:
    using DataType       = CbVector<real>;
    using TransmitterPtr = SPtr<TbTransmitter<DataType>>;

public:
    CreateTransmittersHelper();
    void createTransmitters(const SPtr<Block3D> sblock, const SPtr<Block3D> tblock, int dir, IBlock ib,
                            TransmitterPtr &sender, TransmitterPtr &receiver, std::shared_ptr<vf::parallel::Communicator> comm,
                            TransmitterType tType);

protected:
private:
    std::string generatePoolKey(int srcRank, int srcLevel, int tgtRank, int tgtLevel);
    std::string generateVectorKey(int x1, int x2, int x3, /*int id,*/ int dir, IBlock ib);
    int generateMPITag(int srcLevel, int tgtLevel);
    static unsigned int vKey;
};

#endif

//! \}
