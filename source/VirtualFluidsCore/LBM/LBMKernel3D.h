#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include "LBMSystem.h"
#include "DistributionArray3D.h"
#include "DataSet3D.h"

#include "InterpolationProcessor.h"
#include <MuParser/include/muParser.h>

#include <boost/serialization/serialization.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
class LBMKernel3D;
typedef boost::shared_ptr<LBMKernel3D> LBMKernel3DPtr;

#include "BCProcessor.h"
#include "Block3D.h"

class LBMKernel3D : public boost::enable_shared_from_this<LBMKernel3D>
{
public:
   typedef std::numeric_limits<LBMReal> LBMRealLim;
public:
   LBMKernel3D();
   virtual ~LBMKernel3D();

   virtual LBMKernel3DPtr clone() = 0;

   virtual void calculate() = 0;
   virtual void swapDistributions() = 0;
   virtual double getCallculationTime() = 0;

   virtual void setBCProcessor(BCProcessorPtr bcp);
   virtual BCProcessorPtr getBCProcessor();
   
   virtual void setCollisionFactor(double collFactor);
   virtual double getCollisionFactor() const;
   
   virtual void setGhostLayerWidth(int witdh);
   virtual int  getGhostLayerWidth() const;

   void setDataSet(DataSet3DPtr dataSet);
   DataSet3DPtr getDataSet() const;

   void setForcingX1(LBMReal forcingX1);
   void setForcingX2(LBMReal forcingX2);
   void setForcingX3(LBMReal forcingX3);

   void setForcingX1( const mu::Parser& parser);
   void setForcingX2( const mu::Parser& parser);
   void setForcingX3( const mu::Parser& parser);

   void setForcingX1( const std::string& muParserString);
   void setForcingX2( const std::string& muParserString);
   void setForcingX3( const std::string& muParserString);

   void setIndex(int x1, int x2, int x3);

   void setDeltaT(LBMReal dt);

   bool getCompressible() const;
   void setCompressible(bool val);

   bool getWithForcing() const;
   void setWithForcing(bool val);

   bool getWithSpongeLayer() const;
   void setWithSpongeLayer(bool val);

   void setSpongeLayer(const mu::Parser& parser);
   void setSpongeLayer(const std::string& muParserString);

   void setBlock(Block3DPtr block);
   Block3DPtr getBlock() const;

protected:
   DataSet3DPtr dataSet;
   BCProcessorPtr bcProcessor; 
   LBMReal collFactor;
   int ghostLayerWidth;
   bool compressible;
   
   //forcing 
   bool withForcing;
   mu::Parser muForcingX1;
   mu::Parser muForcingX2;
   mu::Parser muForcingX3;
   int ix1, ix2, ix3;
   LBMReal deltaT;

   //sponge layer
   bool withSpongeLayer;
   mu::Parser muSpongeLayer;

   boost::weak_ptr<Block3D> block;


private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & collFactor;
      ar & ghostLayerWidth;
      ar & compressible;
      ar & withForcing;
      //ar & withSpongeLayer;
      ar & deltaT;
      ar & dataSet;
      ar & bcProcessor;
      ar & ix1;
      ar & ix2;
      ar & ix3;
   }

   void checkFunction(mu::Parser fct);
};

#endif
