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
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Soeren Textor, Sebastian Bindick
//=======================================================================================
#ifndef KDSPLITCANDIDATEMANAGER_H
#define KDSPLITCANDIDATEMANAGER_H

#include <geometry3d/KdTree/KdSplitCandidate.h>

#include <algorithm>
#include <map>
#include <vector>

namespace Kd
{
template <typename T>
class SplitCandidateManager
{
public:
    SplitCandidateManager()

        = default;
    /* ======================================================================================= */
    SplitCandidate<T> &operator[](const std::size_t &i)
    {
#ifdef DEBUG
        return splitCandidatesVec.at(i);
#else
        return splitCandidatesVec[i];
#endif
    }
    /* ======================================================================================= */
    typename std::vector<SplitCandidate<T>>::size_type size() { return splitCandidatesVec.size(); }
    /* ======================================================================================= */
    void add(const T &pos, const int &axis, const int &starting, const int &ending, const int &np)
    {
        typename std::map<T, SplitCandidate<T>>::iterator it = splitCandidates.find(pos);
        if (it != splitCandidates
                      .end()) // split candidate is already available -> increase parameter (starting, ending and np)
        {
            SplitCandidate<T> &sc = it->second;
            sc.np += np;
            sc.starting += starting;
            sc.ending += ending;
        } else // split candidate is not available -> add new split candidate
        {
            this->splitCandidates[pos] = SplitCandidate<T>(axis, pos, starting, ending, np);
        }
    }
    /* ======================================================================================= */
    void createSortedArray()
    {
        splitCandidatesVec.clear();
        typename std::map<T, SplitCandidate<T>>::iterator it;
        for (it = splitCandidates.begin(); it != splitCandidates.end(); ++it)
            splitCandidatesVec.push_back(it->second);
        splitCandidates.clear();
        std::sort(splitCandidatesVec.begin(), splitCandidatesVec.end(), std::less<SplitCandidate<T>>());
    }
    /* ======================================================================================= */

public:
    int objects_starting_outside_left{ 0 };
    int objects_fully_outside_node{ 0 };

private:
    std::map<T, SplitCandidate<T>> splitCandidates;
    std::vector<SplitCandidate<T>> splitCandidatesVec;
};
} // namespace Kd

#endif // KDSPLITCANDIDATEMANAGER_H

//! \}
