#include "EmergencyExitSimulationObserver.h"
#include <mpi/Communicator.h>
#include "Grid3D.h"
#include "MPIIORestartSimulationObserver.h"
#include "UbLogger.h"
#include "UbScheduler.h"
#include <basics/utilities/UbFileInputASCII.h>
#include <basics/utilities/UbFileOutputASCII.h>

EmergencyExitSimulationObserver::EmergencyExitSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                   SPtr<MPIIORestartSimulationObserver> rp, std::shared_ptr<vf::mpi::Communicator> comm)
    : SimulationObserver(grid, s), path(path), rp(rp), comm(comm)
{
    this->path = path + "/exit";
    metafile   = this->path + "/stop.txt";
    if (comm->getProcessID() == comm->getRoot()) {
        // checkMetafile();
        writeMetafile(false);
    }
    comm->barrier();
}
//////////////////////////////////////////////////////////////////////////
EmergencyExitSimulationObserver::~EmergencyExitSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "EmergencyExitSimulationObserver::update:" << step);
}

void EmergencyExitSimulationObserver::collectData(real step)
{
    if (readMetafile()) {
        rp->update((int)step);
        if (comm->getProcessID() == comm->getRoot())
            UBLOG(logINFO, "EmergencyExitSimulationObserver save step: " << step);
        comm->barrier();
        exit(EXIT_SUCCESS);
    }
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::writeMetafile(int /*status*/)
{
    UbFileOutputASCII out(metafile);
    out.writeBool(false);
}
//////////////////////////////////////////////////////////////////////////
bool EmergencyExitSimulationObserver::readMetafile()
{
    UbFileInputASCII in(metafile);
    return in.readBool();
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitSimulationObserver::checkMetafile()
{
    std::ifstream file(metafile.c_str());
    if (!file.is_open()) {
        writeMetafile(false);
        return;
    }
    file.close();
}
