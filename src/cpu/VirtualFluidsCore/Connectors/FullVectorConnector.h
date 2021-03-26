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
//! \file FullVectorConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef FullVectorConnector_H
#define FullVectorConnector_H

#include <vector>

#include "Block3D.h"
#include "RemoteBlock3DConnector.h"

// daten werden in einen vector (dieser befindet sich im transmitter) kopiert
// der vector wird via transmitter uebertragen
// transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
// transmitter sein, der von Transmitter abgeleitet ist ;-)
class FullVectorConnector : public RemoteBlock3DConnector
{
public:
    FullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
                               int sendDir);

    void init() override;

    void fillSendVectors() override;
    void distributeReceiveVectors() override;

protected:
    virtual inline void updatePointers() = 0;
    void fillData();
    void distributeData();
    virtual inline void fillData(vector_type &sdata, int &index, int x1, int x2, int x3) = 0;
    virtual inline void distributeData(vector_type &rdata, int &index, int x1, int x2, int x3) = 0;
    
    int maxX1;
    int maxX2;
    int maxX3;

private:

};

//////////////////////////////////////////////////////////////////////////

#endif // FullVectorConnector_H
