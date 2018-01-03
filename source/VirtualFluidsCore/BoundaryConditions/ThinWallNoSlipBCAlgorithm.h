#ifndef ThinWallNoSlipBCAlgorithm_h__
#define ThinWallNoSlipBCAlgorithm_h__

#include "BCAlgorithm.h"

class DistributionArray3D;

class ThinWallNoSlipBCAlgorithm;
typedef std::shared_ptr<ThinWallNoSlipBCAlgorithm> ThinWallNoSlipBCAlgorithmPtr;

class ThinWallNoSlipBCAlgorithm : public BCAlgorithm
{
public:
   ThinWallNoSlipBCAlgorithm();
   virtual ~ThinWallNoSlipBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(std::shared_ptr<DistributionArray3D> distributions);
   void setPass(int pass);
   void applyBC() override;

protected:
   std::shared_ptr<DistributionArray3D> distributionsTemp;
private:
   int pass;
   LBMReal fTemp[D3Q27System::ENDF + 1];

   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<BCAlgorithm>(*this);
   //}
};
#endif // ThinWallNoSlipBCAlgorithm_h__
