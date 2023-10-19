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
//! \file BCStrategyRegister.h
//! \ingroup BoundarConditions
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BCStrategyRegister_H
#define BCStrategyRegister_H

#include <memory>
#include <map>
#include <mutex>

class BC;
class BCStrategy;

class BCStrategyRegister
{
public:
    BCStrategyRegister(const BCStrategyRegister &) = delete;
    BCStrategyRegister &operator=(const BCStrategyRegister &rhs) = delete;
    static std::shared_ptr<BCStrategyRegister> getInstance();

    virtual ~BCStrategyRegister() = default;

    void setBCStrategy(char bcStrategyKey, std::shared_ptr<BCStrategy> bcStrategy);
    std::shared_ptr<BCStrategy> getBCStrategy(char type);

private:
    BCStrategyRegister() = default;
    static std::mutex instantiation_mutex;
    static std::shared_ptr<BCStrategyRegister> instance;
  
    std::map<char, std::shared_ptr<BCStrategy>> bcMap;
};
#endif