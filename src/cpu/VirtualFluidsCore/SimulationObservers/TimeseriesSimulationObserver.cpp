/*
 *  TimeseriesWriterSimulationObserver.h
 *
 *  Created on: 08.05.2013
 *  Author: uphoff
 */

#include "TimeseriesSimulationObserver.h"

#include <fstream>

#include <mpi/Communicator.h>
#include "Grid3D.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"

TimeseriesSimulationObserver::TimeseriesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<IntegrateValuesHelper> h1,
                                             const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm)
    : SimulationObserver(grid, s), h1(h1), path(path), comm(comm)
{
    if (comm->getProcessID() == comm->getRoot()) {
        std::ofstream ostr;
        // fname = path+"/timeseries/timeseries"+UbSystem::toString(grid->getTimeStep())+".csv";
        fname = path + ".csv";
        UBLOG(logINFO, "TimeseriesWriterSimulationObserver::fname:" << fname);
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string file_path = UbSystem::getPathFromString(fname);
            if (file_path.size() > 0) {
                UbSystem::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << "step;rho;vx;vy;vz;volume\n";
        ostr.close();
        UBLOG(logINFO, "TimeseriesWriterSimulationObserver::Constructor:end");
    }
}
//////////////////////////////////////////////////////////////////////////
TimeseriesSimulationObserver::~TimeseriesSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void TimeseriesSimulationObserver::process(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void TimeseriesSimulationObserver::collectData(real step)
{
    h1->calculateMQ();

    UBLOG(logDEBUG3, "TimeseriesWriterSimulationObserver::update:" << step);

    if (comm->getProcessID() == comm->getRoot()) {
        int istep = static_cast<int>(step);
        std::ofstream ostr;
        real cellsVolume = h1->getCellsVolume();

        real rho    = (h1->getRho()) / cellsVolume;
        real vx     = (h1->getVx1()) / cellsVolume;
        real vy     = (h1->getVx2()) / cellsVolume;
        real vz     = (h1->getVx3()) / cellsVolume;
        real volume = cellsVolume;

        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        ostr << istep << ";" << rho << ";" << vx << ";" << vy << ";" << vz << ";" << volume << "\n";
        ostr.close();
    }
}
