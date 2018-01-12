/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef WriteObjectsCoProcessor_H
#define WriteObjectsCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"

class GbSphere3D;
class Communicator;
class Grid3D;
class UbScheduler;

class WriteObjectsCoProcessor : public  CoProcessor
{
public:
    WriteObjectsCoProcessor();
    WriteObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<Communicator> comm);
   ~WriteObjectsCoProcessor() {}
   void process(double step) override;

   void addGbObject(SPtr<GbSphere3D> sphere);

private:
    std::string path;
    SPtr<Communicator> comm;
    std::vector<SPtr<GbSphere3D> > objects;
};
#endif
