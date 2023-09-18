/*
 *  EmergencyExitSimulationObserver.h
 *
 *  Created on: 05.10.2012
 *  Author: K. Kucher
 */

#ifndef EmergencyExitSimulationObserver_H
#define EmergencyExitSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"

class MPIIORestartSimulationObserver;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;

class EmergencyExitSimulationObserver : public SimulationObserver
{
public:
    EmergencyExitSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                             SPtr<MPIIORestartSimulationObserver> rp, std::shared_ptr<vf::parallel::Communicator> comm);
    ~EmergencyExitSimulationObserver() override;

    void update(real step) override;

protected:
    void collectData(real step);
    void writeMetafile(int status);
    bool readMetafile();
    void checkMetafile();

private:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;
    SPtr<MPIIORestartSimulationObserver> rp;
    std::string metafile;
};

#endif /* EmergencyExitSimulationObserver_H */
