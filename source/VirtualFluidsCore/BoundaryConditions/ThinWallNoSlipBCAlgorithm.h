#ifndef ThinWallNoSlipBCAlgorithm_h__
#define ThinWallNoSlipBCAlgorithm_h__

#include "NoSlipBCAlgorithm.h"

class ThinWallNoSlipBCAlgorithm;
typedef boost::shared_ptr<ThinWallNoSlipBCAlgorithm> ThinWallNoSlipBCAlgorithmPtr;

class ThinWallNoSlipBCAlgorithm : public BCAlgorithm
{
public:
   ThinWallNoSlipBCAlgorithm();
   virtual ~ThinWallNoSlipBCAlgorithm();
   BCAlgorithmPtr clone();
   void addDistributions(DistributionArray3DPtr distributions);
   void setPass(int pass);
protected:
   void applyBC();
   DistributionArray3DPtr distributionsTemp;
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
