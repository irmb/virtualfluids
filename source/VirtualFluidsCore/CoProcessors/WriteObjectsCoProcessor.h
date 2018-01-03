/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef WriteObjectsCoProcessor_H
#define WriteObjectsCoProcessor_H

#include <memory>
#include <string>
#include <memory>

#include "CoProcessor.h"

class GbSphere3D;
class Communicator;
class Grid3D;
class UbScheduler;

class WriteObjectsCoProcessor;
typedef std::shared_ptr<WriteObjectsCoProcessor> WriteObjectsCoProcessorPtr;

class WriteObjectsCoProcessor : public  CoProcessor
{
public:
    WriteObjectsCoProcessor();
    WriteObjectsCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, const std::string& path, std::shared_ptr<Communicator> comm);
   ~WriteObjectsCoProcessor() {}
   void process(double step) override;

   void addGbObject(std::shared_ptr<GbSphere3D> sphere);

private:
    std::string path;
    std::shared_ptr<Communicator> comm;

    std::vector<std::shared_ptr<GbSphere3D> > objects;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {

   }
};
#endif
