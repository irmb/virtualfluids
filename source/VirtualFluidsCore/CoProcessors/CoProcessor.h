#ifndef CoProcessor_H
#define CoProcessor_H

#include <memory>

class Grid3D;
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
protected:
   std::shared_ptr<Grid3D> grid;
   std::shared_ptr<UbScheduler> scheduler;
};
#endif

