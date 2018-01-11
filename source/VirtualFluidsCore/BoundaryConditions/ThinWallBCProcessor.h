#ifndef ThinWallBCProcessor_H
#define ThinWallBCProcessor_H

#include <memory>

#include "BCProcessor.h"

class ThinWallBCProcessor;
typedef std::shared_ptr<ThinWallBCProcessor> ThinWallBCProcessorPtr;

class LBMKernel;

class ThinWallBCProcessor : public BCProcessor
{
public:
   ThinWallBCProcessor();
   ThinWallBCProcessor(std::shared_ptr<LBMKernel> kernel);
   ~ThinWallBCProcessor();
   std::shared_ptr<BCProcessor> clone(std::shared_ptr<LBMKernel> kernel);
   void applyPostCollisionBC();
protected:
private:

};

#endif
