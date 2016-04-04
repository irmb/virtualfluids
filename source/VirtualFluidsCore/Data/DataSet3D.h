#ifndef DataSet3D_h
#define DataSet3D_h

#include <boost/serialization/serialization.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include "DistributionArray3D.h"

class DataSet3D;
typedef boost::shared_ptr<DataSet3D> DataSet3DPtr;

typedef CbArray4D<LBMReal,IndexerX4X3X2X1> AverageValuesArray3D;
typedef boost::shared_ptr< AverageValuesArray3D > AverageValuesArray3DPtr;

typedef CbArray4D<LBMReal, IndexerX4X3X2X1> AverageVelocityArray3D;
typedef boost::shared_ptr< AverageValuesArray3D > AverageVelocityArray3DPtr;

typedef CbArray4D<LBMReal, IndexerX4X3X2X1> AverageFluctuationsArray3D;
typedef boost::shared_ptr< AverageFluctuationsArray3D > AverageFluctuationsArray3DPtr;

typedef CbArray4D<LBMReal, IndexerX4X3X2X1> AverageTriplecorrelationsArray3D;
typedef boost::shared_ptr< AverageValuesArray3D > AverageTriplecorrelationsArray3DPtr;

typedef CbArray4D<LBMReal,IndexerX4X3X2X1> ShearStressValuesArray3D;
typedef boost::shared_ptr< ShearStressValuesArray3D > ShearStressValuesArray3DPtr;

typedef CbArray3D<LBMReal, IndexerX3X2X1> ViscosityArray3D;
typedef boost::shared_ptr< ViscosityArray3D > ViscosityArray3DPtr;

class DataSet3D
{
public:
   DistributionArray3DPtr getFdistributions() const;
   void setFdistributions(DistributionArray3DPtr distributions);

   AverageVelocityArray3DPtr getAverageVelocity() const;
   void setAverageVelocity(AverageVelocityArray3DPtr values);

   AverageFluctuationsArray3DPtr getAverageFluctuations() const;
   void setAverageFluctuations(AverageFluctuationsArray3DPtr values);

   AverageTriplecorrelationsArray3DPtr getAverageTriplecorrelations() const;
   void setAverageTriplecorrelations(AverageTriplecorrelationsArray3DPtr values);



   
   AverageValuesArray3DPtr getAverageValues() const;
   void setAverageValues(AverageValuesArray3DPtr values);
   
   ShearStressValuesArray3DPtr getShearStressValues() const;
   void setShearStressValues(ShearStressValuesArray3DPtr values);
protected:
private:
   DistributionArray3DPtr mFdistributions;
   AverageValuesArray3DPtr mAverageValues;

   AverageVelocityArray3DPtr mAverageVelocity;
   AverageFluctuationsArray3DPtr mAverageFluktuations;
   AverageTriplecorrelationsArray3DPtr mAverageTriplecorrelations;

   ShearStressValuesArray3DPtr mShearStressValues;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & mFdistributions;
      ar & mAverageValues;
      ar & mShearStressValues;
      ar & mAverageVelocity;
      ar & mAverageFluktuations;
      ar & mAverageTriplecorrelations;
   }
};

inline DistributionArray3DPtr DataSet3D::getFdistributions() const
{ 
   return mFdistributions; 
}

inline void DataSet3D::setFdistributions(DistributionArray3DPtr distributions)
{ 
   mFdistributions = distributions; 
}

inline AverageValuesArray3DPtr DataSet3D::getAverageValues() const
{ 
   return mAverageValues; 
}

inline void DataSet3D::setAverageValues(AverageValuesArray3DPtr values)
{ 
   mAverageValues = values; 
}

inline AverageVelocityArray3DPtr DataSet3D::getAverageVelocity() const
{
   return mAverageVelocity;
}

inline void DataSet3D::setAverageVelocity(AverageVelocityArray3DPtr values)
{
   mAverageVelocity = values;
}

inline AverageFluctuationsArray3DPtr DataSet3D::getAverageFluctuations() const
{
   return mAverageFluktuations;
}

inline void DataSet3D::setAverageFluctuations(AverageFluctuationsArray3DPtr values)
{
   mAverageFluktuations = values;
}

inline AverageTriplecorrelationsArray3DPtr DataSet3D::getAverageTriplecorrelations() const
{
   return mAverageTriplecorrelations;
}

inline void DataSet3D::setAverageTriplecorrelations(AverageTriplecorrelationsArray3DPtr values)
{
   mAverageTriplecorrelations = values;
}

inline ShearStressValuesArray3DPtr DataSet3D::getShearStressValues() const
{ 
   return mShearStressValues; 
}

inline void DataSet3D::setShearStressValues(ShearStressValuesArray3DPtr values)
{ 
   mShearStressValues = values; 
}
#endif
