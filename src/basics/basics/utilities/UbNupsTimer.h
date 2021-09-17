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
//! \file UbNupsTimer.h
//! \ingroup utilities
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBNUPSTIMER_H
#define UBNUPSTIMER_H

#include <basics/utilities/UbTiming.h>
#include <sstream>
#include <vector>

/*=========================================================================*/
/*  UbNupsTimer                                                             */
/*                                                                         */
/**
This Class provides the base for ...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 01.11.04
*/
class UbNupsTimer : public UbTiming
{
public:
    UbNupsTimer() : UbTiming()
    {
        mTempNodes = 0.0;
        mNofNodes.resize(0);
        mDurations.resize(0);
    }
    /*==========================================================*/
    UbNupsTimer(std::string name) : UbTiming(name)
    {
        mNofNodes.resize(0);
        mDurations.resize(0);
        mTempNodes = 0.0;
    }
    /*==========================================================*/
    void initTiming() override
    {
        UbTiming::initTiming();
        mNofNodes.resize(0);
        mDurations.resize(0);
        mTempNodes = 0.0;
    }
    /*==========================================================*/
    void startNUPSTiming(double nofNodes)
    {
        mTempNodes = nofNodes;
        UbTiming::startTiming();
    }
    /*==========================================================*/
    void endTiming() override
    {
        UbTiming::endTiming();
        // save #node and time informations
        mNofNodes.push_back(mTempNodes);
        mDurations.push_back(UbTiming::getDuration());
        // reset internal timecounter
        UbTiming::initTiming();
    }
    /*==========================================================*/
    double getAverageNups()
    {
        double averageNups = 0.0;
        for (int i = 0; i < (int)mNofNodes.size(); i++)
            averageNups += mNofNodes.at(i) / mDurations.at(i);

        return averageNups / (double)mNofNodes.size();
    }
    /*==========================================================*/
    double getSumOfDuration()
    {
        double duration = 0.0;
        for (int i = 0; i < (int)mDurations.size(); i++)
            duration += mDurations.at(i);
        return duration;
    }
    /*==========================================================*/
    std::string getNupsString()
    {
        std::stringstream ss;
        ss << "saved nups informations" << std::endl;
        for (int i = 0; i < (int)mNofNodes.size(); i++)
            ss << mNofNodes.at(i) << "nodes/" << mDurations.at(i) << "sec=" << mNofNodes.at(i) / mDurations.at(i)
               << "nups\n";
        return ss.str();
    }

protected:
private:
    std::vector<double> mNofNodes;
    std::vector<double> mDurations;

    double mTempNodes;
};

#endif
