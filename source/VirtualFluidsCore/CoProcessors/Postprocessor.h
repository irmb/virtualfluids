#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>

#include <Grid3D.h>
#include <basics/utilities/UbScheduler.h>

class Postprocessor;
typedef boost::shared_ptr<Postprocessor> PostprocessorPtr;

class Postprocessor
{
public:
   Postprocessor()
   {

   }

   Postprocessor(Grid3DPtr grid, UbSchedulerPtr s)
      : grid(grid)
      , scheduler(s)
   {
      connection = grid->connect(boost::bind(&Postprocessor::update, this, _1));
   }

   virtual ~Postprocessor()
   {
      grid->disconnect(connection);
   }

   virtual void update(double step) = 0;

   virtual void disconnect()
   {
      grid->disconnect(connection);
   }

   virtual void reconnect(Grid3DPtr grid)
   {
      this->grid = grid;
      this->grid->disconnect(connection);
      connection = this->grid->connect(boost::bind(&Postprocessor::update, this, _1));
   }

protected:
   Grid3DPtr grid;
   UbSchedulerPtr scheduler;
private:
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

