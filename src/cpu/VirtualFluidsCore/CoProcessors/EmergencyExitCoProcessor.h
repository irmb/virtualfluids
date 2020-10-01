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
class Communicator;
class Grid3D;
class UbScheduler;

class EmergencyExitCoProcessor : public CoProcessor
{
public:
    EmergencyExitCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<MPIIORestartCoProcessor> rp, SPtr<Communicator> comm);
    virtual ~EmergencyExitCoProcessor();

    void process(double step) override;

protected:
    void collectData(double step);
    void writeMetafile(int status);
    bool readMetafile();
    void checkMetafile();

private:
    std::string path;
    SPtr<Communicator> comm;
    SPtr<MPIIORestartCoProcessor> rp;
    std::string metafile;
};


#endif /* EmergencyExitCoProcessor_H */
