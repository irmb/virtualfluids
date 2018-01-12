#ifndef CoProcessor_H
#define CoProcessor_H

#include <PointerDefinitions.h>

class Grid3D;
class UbScheduler;

class CoProcessor
{
public:
   CoProcessor();
   CoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
   virtual ~CoProcessor();
   virtual void process(double step) = 0;
protected:
   SPtr<Grid3D> grid;
   SPtr<UbScheduler> scheduler;
};
#endif

