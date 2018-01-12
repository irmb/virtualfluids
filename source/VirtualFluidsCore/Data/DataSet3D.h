#ifndef DataSet3D_h
#define DataSet3D_h

#include <PointerDefinitions.h>

#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
#include "DistributionArray3D.h"

typedef CbArray4D<LBMReal,IndexerX4X3X2X1> AverageValuesArray3D;
typedef CbArray4D<LBMReal,IndexerX4X3X2X1> ShearStressValuesArray3D;
typedef CbArray3D<LBMReal, IndexerX3X2X1> RelaxationFactorArray3D;

class DataSet3D
{
public:
   SPtr<DistributionArray3D> getFdistributions() const;
   void setFdistributions(SPtr<DistributionArray3D> distributions);

   SPtr<AverageValuesArray3D> getAverageDencity() const;
   void setAverageDencity(SPtr<AverageValuesArray3D> values);

   SPtr<AverageValuesArray3D> getAverageVelocity() const;
   void setAverageVelocity(SPtr<AverageValuesArray3D> values);

   SPtr<AverageValuesArray3D> getAverageFluctuations() const;
   void setAverageFluctuations(SPtr<AverageValuesArray3D> values);

   SPtr<AverageValuesArray3D> getAverageTriplecorrelations() const;
   void setAverageTriplecorrelations(SPtr<AverageValuesArray3D> values);
   
   SPtr<AverageValuesArray3D> getAverageValues() const;
   void setAverageValues(SPtr<AverageValuesArray3D> values);
   
   SPtr<ShearStressValuesArray3D> getShearStressValues() const;
   void setShearStressValues(SPtr<ShearStressValuesArray3D> values);

   SPtr<RelaxationFactorArray3D> getRelaxationFactor() const;
   void setRelaxationFactor(SPtr<RelaxationFactorArray3D> values);
protected:
private:
   SPtr<DistributionArray3D> fdistributions;
   SPtr<AverageValuesArray3D> averageValues;

   SPtr<AverageValuesArray3D> averageDencity;
   SPtr<AverageValuesArray3D> averageVelocity;
   SPtr<AverageValuesArray3D> averageFluktuations;
   SPtr<AverageValuesArray3D> averageTriplecorrelations;

   SPtr<ShearStressValuesArray3D> shearStressValues;

   SPtr<RelaxationFactorArray3D> relaxationFactor;

};

inline SPtr<DistributionArray3D> DataSet3D::getFdistributions() const
{ 
   return fdistributions; 
}

inline void DataSet3D::setFdistributions(SPtr<DistributionArray3D> distributions)
{ 
   fdistributions = distributions; 
}

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageValues() const
{ 
   return averageValues; 
}

inline void DataSet3D::setAverageValues(SPtr<AverageValuesArray3D> values)
{ 
   averageValues = values; 
}

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageDencity() const
{
   return averageDencity;
}

inline void DataSet3D::setAverageDencity(SPtr<AverageValuesArray3D> values)
{
   averageDencity = values;
}

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageVelocity() const
{
   return averageVelocity;
}

inline void DataSet3D::setAverageVelocity(SPtr<AverageValuesArray3D> values)
{
   averageVelocity = values;
}

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageFluctuations() const
{
   return averageFluktuations;
}

inline void DataSet3D::setAverageFluctuations(SPtr<AverageValuesArray3D> values)
{
   averageFluktuations = values;
}

inline SPtr<AverageValuesArray3D> DataSet3D::getAverageTriplecorrelations() const
{
   return averageTriplecorrelations;
}

inline void DataSet3D::setAverageTriplecorrelations(SPtr<AverageValuesArray3D> values)
{
   averageTriplecorrelations = values;
}

inline SPtr<ShearStressValuesArray3D> DataSet3D::getShearStressValues() const
{ 
   return shearStressValues; 
}

inline void DataSet3D::setShearStressValues(SPtr<ShearStressValuesArray3D> values)
{ 
   shearStressValues = values; 
}

inline SPtr<RelaxationFactorArray3D> DataSet3D::getRelaxationFactor() const
{
   return relaxationFactor;
}

inline void DataSet3D::setRelaxationFactor(SPtr<RelaxationFactorArray3D> values)
{
   relaxationFactor = values;
}
#endif
