#ifndef BC_PROCESSSOR_H
#define BC_PROCESSSOR_H

#include <memory>
#include <vector>

#include <boost/serialization/base_object.hpp>

class BCProcessor;
typedef std::shared_ptr<BCProcessor> BCProcessorPtr;

class BCArray3D;
class BCAlgorithm;
class ILBMKernel;

class BCProcessor
{
public:
   BCProcessor();
   BCProcessor(std::shared_ptr<ILBMKernel> kernel);
   virtual ~BCProcessor();
   virtual std::shared_ptr<BCArray3D> getBCArray();
   virtual void setBCArray(std::shared_ptr<BCArray3D> bcarray);
   virtual BCProcessorPtr clone(std::shared_ptr<ILBMKernel> kernel);

   void addBC(std::shared_ptr<BCAlgorithm> bc);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
   void clearBC();
protected:
   std::vector<std::shared_ptr<BCAlgorithm> > preBC;
   std::vector<std::shared_ptr<BCAlgorithm> > postBC;
   std::shared_ptr<BCArray3D> bcArray;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & bcArray;
      //ar & preBC;
      //ar & postBC;
   }
};

#endif
