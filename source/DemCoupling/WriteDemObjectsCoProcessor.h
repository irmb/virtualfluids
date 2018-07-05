/*
*  Author: K. Kutscher
*  mail: kutscher@irmb.tu-bs.de
*/
#ifndef WriteDemObjectsCoProcessor_H
#define WriteDemObjectsCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class DemCoProcessor;
class WbWriter;

class WriteDemObjectsCoProcessor : public  CoProcessor
{
public:
    WriteDemObjectsCoProcessor();
    WriteDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<DemCoProcessor> demCoProcessor, SPtr<Communicator> comm);
   ~WriteDemObjectsCoProcessor() {}
   void process(double step) override;

private:
    std::string path;
    WbWriter* writer;
    SPtr<Communicator> comm;
    SPtr<DemCoProcessor> demCoProcessor;
};
#endif
