#ifndef ThinWallBCProcessor_H
#define ThinWallBCProcessor_H

#include <PointerDefinitions.h>

#include "BCProcessor.h"

class LBMKernel;

class ThinWallBCProcessor : public BCProcessor
{
public:
   ThinWallBCProcessor();
   ThinWallBCProcessor(SPtr<LBMKernel> kernel);
   ~ThinWallBCProcessor();
   SPtr<BCProcessor> clone(SPtr<LBMKernel> kernel);
   void applyPostCollisionBC();
protected:
private:

};

#endif
