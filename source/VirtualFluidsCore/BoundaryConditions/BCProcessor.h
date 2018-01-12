#ifndef BC_PROCESSSOR_H
#define BC_PROCESSSOR_H

#include <PointerDefinitions.h>
#include <vector>

class BCArray3D;
class BCAlgorithm;
class ILBMKernel;

class BCProcessor
{
public:
   BCProcessor();
   BCProcessor(SPtr<ILBMKernel> kernel);
   virtual ~BCProcessor();
   virtual SPtr<BCArray3D> getBCArray();
   virtual void setBCArray(SPtr<BCArray3D> bcarray);
   virtual SPtr<BCProcessor> clone(SPtr<ILBMKernel> kernel);

   void addBC(SPtr<BCAlgorithm> bc);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
   void clearBC();
protected:
   std::vector<SPtr<BCAlgorithm> > preBC;
   std::vector<SPtr<BCAlgorithm> > postBC;
   SPtr<BCArray3D> bcArray;

private:

};

#endif
