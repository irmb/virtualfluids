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
//! \file MbSmartPtrBase.cpp
//! \ingroup memory
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#include <basics/memory/MbSmartPtrBase.h>

using namespace std;

bool MbSmartPtrBase::addRef(void *ptr)
{
    MbSmartPtrBaseMap::getInstance()->getMap()[ptr]++;
    return true;
}
//-------------------------------------------------
bool MbSmartPtrBase::releaseRef(void *ptr)
{
    map<void *, int> &ptrMap       = MbSmartPtrBaseMap::getInstance()->getMap();
    map<void *, int>::iterator pos = ptrMap.find(ptr);

    if (pos != ptrMap.end()) {
        pos->second--;

        if (pos->second == 0) {
            ptrMap.erase(pos);
            return true;
        }
    }
    return false;
}
//-------------------------------------------------
bool MbSmartPtrBase::removeFromGC(void *ptr) const
{
    if (MbSmartPtrBaseMap::getInstance()->getMap().erase(ptr))
        return true;
    return false;
}
//-------------------------------------------------
int MbSmartPtrBase::ref_count(void *ptr) const
{
    map<void *, int> &ptrMap       = MbSmartPtrBaseMap::getInstance()->getMap();
    map<void *, int>::iterator pos = ptrMap.find(ptr);

    if (pos != ptrMap.end())
        return pos->second;
    else
        return 0;
}
