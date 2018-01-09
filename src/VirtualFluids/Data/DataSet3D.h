#ifndef DataSet3D_h
#define DataSet3D_h

#include <boost/serialization/serialization.hpp>
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include "DistributionArray3D.h"

class DataSet3D;
typedef std::shared_ptr<DataSet3D> DataSet3DPtr;

typedef CbArray4D<LBMReal,IndexerX4X3X2X1> AverageValuesArray3D;
typedef std::shared_ptr< AverageValuesArray3D > AverageValuesArray3DPtr;

typedef CbArray4D<LBMReal,IndexerX4X3X2X1> ShearStressValuesArray3D;
typedef std::shared_ptr< ShearStressValuesArray3D > ShearStressValuesArray3DPtr;

typedef CbArray3D<LBMReal, IndexerX3X2X1> RelaxationFactorArray3D;
typedef std::shared_ptr< RelaxationFactorArray3D > RelaxationFactorArray3DPtr;

class DataSet3D
{
public:
   DistributionArray3DPtr getFdistributions() const;
   void setFdistributions(DistributionArray3DPtr distributions);

   AverageValuesArray3DPtr getAverageDencity() const;
   void setAverageDencity(AverageValuesArray3DPtr values);

   AverageValuesArray3DPtr getAverageVelocity() const;
   void setAverageVelocity(AverageValuesArray3DPtr values);

   AverageValuesArray3DPtr getAverageFluctuations() const;
   void setAverageFluctuations(AverageValuesArray3DPtr values);

   AverageValuesArray3DPtr getAverageTriplecorrelations() const;
   void setAverageTriplecorrelations(AverageValuesArray3DPtr values);
   
   AverageValuesArray3DPtr getAverageValues() const;
   void setAverageValues(AverageValuesArray3DPtr values);
   
   ShearStressValuesArray3DPtr getShearStressValues() const;
   void setShearStressValues(ShearStressValuesArray3DPtr values);

   RelaxationFactorArray3DPtr getRelaxationFactor() const;
   void setRelaxationFactor(RelaxationFactorArray3DPtr values);
protected:
private:
   DistributionArray3DPtr fdistributions;
   AverageValuesArray3DPtr averageValues;

   AverageValuesArray3DPtr averageDencity;
   AverageValuesArray3DPtr averageVelocity;
   AverageValuesArray3DPtr averageFluktuations;
   AverageValuesArray3DPtr averageTriplecorrelations;

   ShearStressValuesArray3DPtr shearStressValues;

   RelaxationFactorArray3DPtr relaxationFactor;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & fdistributions;
      ar & averageValues;
      ar & shearStressValues;
      ar & averageDencity;
      ar & averageVelocity;
      ar & averageFluktuations;
      ar & averageTriplecorrelations;
      ar & relaxationFactor;
   }
};

inline DistributionArray3DPtr DataSet3D::getFdistributions() const
{ 
   return fdistributions; 
}

inline void DataSet3D::setFdistributions(DistributionArray3DPtr distributions)
{ 
   fdistributions = distributions; 
}

inline AverageValuesArray3DPtr DataSet3D::getAverageValues() const
{ 
   return averageValues; 
}

inline void DataSet3D::setAverageValues(AverageValuesArray3DPtr values)
{ 
   averageValues = values; 
}

inline AverageValuesArray3DPtr DataSet3D::getAverageDencity() const
{
   return averageDencity;
}

inline void DataSet3D::setAverageDencity(AverageValuesArray3DPtr values)
{
   averageDencity = values;
}

inline AverageValuesArray3DPtr DataSet3D::getAverageVelocity() const
{
   return averageVelocity;
}

inline void DataSet3D::setAverageVelocity(AverageValuesArray3DPtr values)
{
   averageVelocity = values;
}

inline AverageValuesArray3DPtr DataSet3D::getAverageFluctuations() const
{
   return averageFluktuations;
}

inline void DataSet3D::setAverageFluctuations(AverageValuesArray3DPtr values)
{
   averageFluktuations = values;
}

inline AverageValuesArray3DPtr DataSet3D::getAverageTriplecorrelations() const
{
   return averageTriplecorrelations;
}

inline void DataSet3D::setAverageTriplecorrelations(AverageValuesArray3DPtr values)
{
   averageTriplecorrelations = values;
}

inline ShearStressValuesArray3DPtr DataSet3D::getShearStressValues() const
{ 
   return shearStressValues; 
}

inline void DataSet3D::setShearStressValues(ShearStressValuesArray3DPtr values)
{ 
   shearStressValues = values; 
}

inline RelaxationFactorArray3DPtr DataSet3D::getRelaxationFactor() const
{
   return relaxationFactor;
}

inline void DataSet3D::setRelaxationFactor(RelaxationFactorArray3DPtr values)
{
   relaxationFactor = values;
}
#endif
