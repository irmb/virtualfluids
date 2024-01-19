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
//! \addtogroup utilities
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller, Jan Hegewald
//=======================================================================================
#ifndef UBSCHEDULER_H
#define UBSCHEDULER_H

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <basics/utilities/UbComparators.h>
#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbMath.h>
#include <basics/utilities/UbSystem.h>

//////////////////////////////////////////////////////////////////////////
//!
//! \brief A class implements scheduling.
//! \details This class is not thread save.
//!
//////////////////////////////////////////////////////////////////////////

class UbScheduler
{
public:
    class UbSchedule
    {
        friend class UbScheduler;

    public:
        UbSchedule() : step(Ub::inf), begin(Ub::inf), end(Ub::inf) {}
        UbSchedule(const double &step, const double &begin = 0.0, const double &end = Ub::inf)
            : step(step), begin(begin), end(end)
        {
        }
        double getStep() const { return this->step; }
        double getBegin() const { return this->begin; }
        double getEnd() const { return this->end; }

        /*==========================================================*/
        std::string toString()
        {
            std::stringstream text;
            text << *this;
            return text.str();
        }
        /*==========================================================*/
        friend inline std::ostream &operator<<(std::ostream &os, const UbSchedule &schedule)
        {
            os << "Schedule[start,end,step]=[" << schedule.begin << ", " << schedule.end << ", " << schedule.step
               << "]";
            return os;
        }

    private:
        double step, begin, end;
    };

public:
    UbScheduler() { this->initVals(); }
    /*==========================================================*/
    UbScheduler(const double &step, const double &begin = 0, const double &end = Ub::inf)
    {
        this->initVals();
        this->addSchedule(step, begin, end);
    }
    /*==========================================================*/
    UbScheduler(const UbSchedule &schedule)
    {
        this->initVals();
        this->addSchedule(schedule);
    }
    /*==========================================================*/
    virtual ~UbScheduler() = default;
    /*==========================================================*/
    inline void addSchedule(const UbSchedule &schedule)
    {
        this->addSchedule(schedule.step, schedule.begin, schedule.end);
    }
    /*==========================================================*/
    bool addSchedule(const double &step, const double &begin, double end)
    {
        if (UbMath::zero(step) || begin > end) {
            std::cerr << "UbScheduler::addSchedule - invalid Schedule:\n\t" << UbSchedule(step, begin, end)
                      << std::endl;
            return false;
        }

        if (UbMath::less(end, (double)Ub::inf)) {
            // es kann vorkommen, dass man mit dem intervall nicht genau auf den letzten wert kommt
            //(z.B. step=2; start=0; end=9; -> ende wird angepasst)
            // also wenn end-begin>Ub::inf ist, dann geht es halt nicht.. ein cast in long double half hier nichts
            double multiplier = 0.0;
            double fractpart  = modf((end - begin) / step, &multiplier);
            if (!UbMath::zero(fractpart)) {
                // tmp-speicherung (fuer cerr)
                fractpart = end;
                // neues ende
                end = begin + multiplier * step;

                std::cerr << "Warning: UbScheduler::addSchedule - "
                          << "end of schedule was adapted to intervall \n\t"
                          << "from " << UbSchedule(step, begin, fractpart) << " to " << UbSchedule(step, begin, end)
                          << std::endl;
            }
        }

        // nu aber:
        schedules.emplace_back(step, begin, end);

        if (end > maxT)
            maxT = end;

        double potentialDueTime;
        if (calcNextDueTimeForSchedule(schedules.back(), lastUsedT, potentialDueTime) &&
            potentialDueTime < nextDueTime) {
            nextDueTime = potentialDueTime;
        }

        return true;
    }
    /*==========================================================*/
    // returns true if scheduler contains schedules
    bool hasSchedules() const { return !schedules.empty(); }
    /*==========================================================*/
    // time bei dem das letzte mal isDue(time) true war
    double getLastDueTime() const { return lastDueTime; }
    /*==========================================================*/
    // time bei dem das naechste mal isDue(time) true ergibt
    double getNextDueTime() const { return nextDueTime; }
    /*==========================================================*/
    // maxDueTime (maxTime der Schedules!
    double getMaxDueTime() const { return this->maxT; }
    /*==========================================================*/
    bool isDue(const double &t)
    {
        lastUsedT = t;
        if (UbMath::greaterEqual(t, nextDueTime)) {
            // groesser maxT is nicht
            if (UbMath::greater(t, maxT))
                return false;

            // temp var
            double actDueTime = nextDueTime;

            // um Suche nach nextDueTime bei "Zukunfts-t" zu optimieren, setzt man die "start"-suchzeit auf "t-1":
            nextDueTime = t - 1; // t-1 deshlab, damit falls z.B. while Schleife nicht durchlaufen wird
                                 // die folgende if Abfrage nicht faelschlicher Weise true ist!
            while (UbMath::greaterEqual(t, nextDueTime) && !UbMath::equal(nextDueTime, maxT)) {
                double tmpNextDueTime = maxT, potentialDueTime = -1.0;
                for (std::size_t i = 0; i < schedules.size(); i++) {
                    if (calcNextDueTimeForSchedule(schedules[i], nextDueTime, potentialDueTime) &&
                        potentialDueTime < tmpNextDueTime) {
                        assert(nextDueTime < potentialDueTime);
                        tmpNextDueTime = potentialDueTime;
                    }
                }
                actDueTime  = nextDueTime;
                nextDueTime = tmpNextDueTime;
            }

            // wenn t = der aktuuellen oder gar schon der naechstmoeglichen ist (hierbei wurde
            // zuvor actDueTime und nextDueTime ggf. angepasst)
            // Bsp.: nextDuTime war 5, aber fuer t=400 gilt andere schedule -> Bsp actDue=350 und nextDue 405
            if (UbMath::equal(t, actDueTime) || UbMath::equal(t, nextDueTime)) {
                lastDueTime = t;
                return true;
            }
        } else if (UbMath::lessEqual(t, lastDueTime)) {
            if (UbMath::equal(t, lastDueTime))
                return true; // braucht man, wenn man fuer dasselbe t isDue(t) aufruft
            else {
                // Fall: Zeit liegt faktisch in der Vergangenheit -> neu initialsisieren
                double tmpNextDueTime = maxT, potentialDueTime = -1.0;
                for (size_t i = 0; i < schedules.size(); i++) {
                    if (calcNextDueTimeForSchedule(schedules[i], t - 1, potentialDueTime) &&
                        potentialDueTime < tmpNextDueTime) {
                        tmpNextDueTime = potentialDueTime;
                    }
                }
                nextDueTime = tmpNextDueTime;

                return UbMath::equal(t, nextDueTime);
            }
        }

        return false;
    }
    /*==========================================================*/
    inline double getMinBegin() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::min_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getBegin))
            ->getBegin();
    }
    /*==========================================================*/
    inline double getMaxBegin() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::max_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getBegin))
            ->getBegin();
    }
    /*==========================================================*/
    inline double getMinEnd() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::min_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getEnd))
            ->getEnd();
    }
    /*==========================================================*/
    inline double getMaxEnd() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::max_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getEnd))
            ->getEnd();
    }
    /*==========================================================*/
    inline double getMinStep() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::min_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getStep))
            ->getStep();
    }
    /*==========================================================*/
    inline double getMaxStep() const
    {
        if (schedules.empty())
            return Ub::inf;
        return std::max_element(schedules.begin(), schedules.end(), UbComparators::membercomp(&UbSchedule::getStep))
            ->getStep();
    }
    /*==========================================================*/
    inline std::string toString() const
    {
        std::stringstream text;
        text << *this;
        return text.str();
    }
    /*==========================================================*/
    friend inline std::ostream &operator<<(std::ostream &os, const UbScheduler &scheduler)
    {
        os << "UbScheduler\n";
        os << "Schedule |       start       |        end        |     intervall     " << std::endl;
        for (std::size_t i = 0; i < scheduler.schedules.size(); i++)
            os << std::setw(9) << i << "|" << std::setw(19) << scheduler.schedules[i].getBegin() << "|" << std::setw(19)
               << scheduler.schedules[i].getEnd() << "|" << std::setw(19) << scheduler.schedules[i].getStep()
               << std::endl;
        return os;
    }

protected:
    /*==========================================================*/
    void initVals()
    {
        lastUsedT   = -Ub::inf;
        lastDueTime = -Ub::inf;
        nextDueTime = Ub::inf;
        maxT        = -Ub::inf;
    }
    /*==========================================================*/
    // calculates next due time for a schedule
    // with  nextDueTime > searchStart
    bool calcNextDueTimeForSchedule(const UbSchedule &schedule, const double &searchStart, double &nextDueTime)
    {
        if (UbMath::greater(searchStart, schedule.end))
            return false;
        else if (UbMath::less(searchStart, schedule.begin))
            nextDueTime = schedule.begin;
        else {
            nextDueTime = schedule.begin + ((int)((searchStart - schedule.begin) / schedule.step) + 1) * schedule.step;
            if (UbMath::less(nextDueTime, searchStart) || UbMath::greater(nextDueTime, schedule.end)) {
                return false;
            }
        }
        return true;
    }

protected:
    double lastUsedT;
    double lastDueTime;
    double nextDueTime;
    double maxT;

    std::vector<UbSchedule> schedules;
};

using UbSchedule = UbScheduler::UbSchedule;

#endif // UBSCHEDULER_H

// int main(int argc, char** argv)
//{
//    UbScheduler writeSchedule;
////    writeSchedule.addSchedule(0,2000,100);
////    writeSchedule.addSchedule(3005,4500,300);
////    writeSchedule.addSchedule(0,10,1);
////    writeSchedule.addSchedule(0,100001,100);
//    writeSchedule.addSchedule(0,2,1);
//    writeSchedule.addSchedule(0,100001,200);
//
//    for(int t = 0; t < 1001; t++)
//    {
//        if(writeSchedule.isDue(t))
//        {
//            cout<<"due@ "<<t<<endl;
//        }
//    }
//    return 0;
//}

//! \}
