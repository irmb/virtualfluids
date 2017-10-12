#ifndef CoProcessor_H
#define CoProcessor_H

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>

#include <Grid/Grid3D.h>
#include <basics/utilities/UbScheduler.h>

class CoProcessor;
typedef boost::shared_ptr<CoProcessor> CoProcessorPtr;

class CoProcessor
{
public:
   CoProcessor()
   {

   }

   CoProcessor(Grid3DPtr grid, UbSchedulerPtr s)
      : grid(grid)
      , scheduler(s)
   {
      connection = grid->connect(boost::bind(&CoProcessor::process, this, _1));
   }

   virtual ~CoProcessor()
   {
      grid->disconnect(connection);
   }

   virtual void process(double step) = 0;

   virtual void disconnect()
   {
      grid->disconnect(connection);
   }

   virtual void reconnect(Grid3DPtr grid)
   {
      this->grid = grid;
      this->grid->disconnect(connection);
      connection = this->grid->connect(boost::bind(&CoProcessor::process, this, _1));
   }

protected:
   Grid3DPtr grid;
   UbSchedulerPtr scheduler;
   Grid3D::connection_t connection;
private:
   

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

