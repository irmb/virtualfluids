/*
 *  EmergencyExitCoProcessor.h
 *
 *  Created on: 05.10.2012
 *  Author: K. Kucher
 */

#ifndef EmergencyExitCoProcessor_H
#define EmergencyExitCoProcessor_H

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

class MPIIORestartCoProcessor;
namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;

class EmergencyExitCoProcessor : public CoProcessor
{
public:
    EmergencyExitCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                             SPtr<MPIIORestartCoProcessor> rp, std::shared_ptr<vf::mpi::Communicator> comm);
    ~EmergencyExitCoProcessor() override;

    void process(real step) override;

protected:
    void collectData(real step);
    void writeMetafile(int status);
    bool readMetafile();
    void checkMetafile();

private:
    std::string path;
    std::shared_ptr<vf::mpi::Communicator> comm;
    SPtr<MPIIORestartCoProcessor> rp;
    std::string metafile;
};

#endif /* EmergencyExitCoProcessor_H */
