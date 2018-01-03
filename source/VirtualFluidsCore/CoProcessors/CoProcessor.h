#ifndef CoProcessor_H
#define CoProcessor_H

#include <memory>

#include <boost/serialization/serialization.hpp>

#include "Grid3D.h"

class UbScheduler;

class CoProcessor;
typedef std::shared_ptr<CoProcessor> CoProcessorPtr;

class CoProcessor
{
public:
    CoProcessor();
    CoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s);
    virtual ~CoProcessor();

   virtual void process(double step) = 0;

   virtual void disconnect();
   virtual void reconnect(std::shared_ptr<Grid3D> grid);

protected:
   std::shared_ptr<Grid3D> grid;
   std::shared_ptr<UbScheduler> scheduler;
   Grid3D::connection_t connection;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      //ar & grid;
      ar & scheduler;
   }
};
#endif

