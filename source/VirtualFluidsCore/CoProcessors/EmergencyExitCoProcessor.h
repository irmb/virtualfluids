/*
 *  EmergencyExitCoProcessor.h
 *
 *  Created on: 05.10.2012
 *  Author: K. Kucher
 */

#ifndef EmergencyExitCoProcessor_H
#define EmergencyExitCoProcessor_H

#include <memory>
#include <string>

#include "CoProcessor.h"

class MPIIORestartCoProcessor;
class Communicator;
class Grid3D;
class UbScheduler;

class EmergencyExitCoProcessor;
typedef std::shared_ptr<EmergencyExitCoProcessor> EmergencyExitCoProcessorPtr;

class EmergencyExitCoProcessor : public CoProcessor
{
public:
    EmergencyExitCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, const std::string& path, std::shared_ptr<MPIIORestartCoProcessor> rp, std::shared_ptr<Communicator> comm);
    virtual ~EmergencyExitCoProcessor();

    void process(double step) override;

protected:
    void collectData(double step);
    void writeMetafile(int status);
    bool readMetafile();
    void checkMetafile();

private:
    std::string path;
    std::shared_ptr<Communicator> comm;
    std::shared_ptr<MPIIORestartCoProcessor> rp;
    std::string metafile;
};


#endif /* EmergencyExitCoProcessor_H */
