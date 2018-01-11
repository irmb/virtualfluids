#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include <memory>

#include "LBMSystem.h"

#include <MuParser/include/muParser.h>

#include "ILBMKernel.h"

class LBMKernel;
typedef std::shared_ptr<LBMKernel> LBMKernelPtr;

class BCProcessor;
class DataSet3D;
class Block3D;

class LBMKernel : public ILBMKernel, public std::enable_shared_from_this<LBMKernel>
{
public:
    typedef std::numeric_limits<LBMReal> LBMRealLim;
public:
    LBMKernel();
    virtual ~LBMKernel();

    virtual LBMKernelPtr clone() = 0;

    virtual void calculate() = 0;
    virtual double getCalculationTime() = 0;

    void setBCProcessor(std::shared_ptr<BCProcessor> bcp);
    std::shared_ptr<BCProcessor> getBCProcessor() const;

    void setCollisionFactor(double collFactor);
    double getCollisionFactor() const;

    void setGhostLayerWidth(int witdh);
    int  getGhostLayerWidth() const;

    void setDataSet(std::shared_ptr<DataSet3D> dataSet);
    std::shared_ptr<DataSet3D> getDataSet() const;

    void setForcingX1(LBMReal forcingX1);
    void setForcingX2(LBMReal forcingX2);
    void setForcingX3(LBMReal forcingX3);

    void setForcingX1(const mu::Parser& parser);
    void setForcingX2(const mu::Parser& parser);
    void setForcingX3(const mu::Parser& parser);

    void setForcingX1(const std::string& muParserString);
    void setForcingX2(const std::string& muParserString);
    void setForcingX3(const std::string& muParserString);

    void setIndex(int x1, int x2, int x3);

    LBMReal getDeltaT() const;
    void setDeltaT(LBMReal dt);

    bool getCompressible() const;
    void setCompressible(bool val);

    bool getWithForcing() const;
    void setWithForcing(bool val);

    bool getWithSpongeLayer() const;
    void setWithSpongeLayer(bool val);

    void setSpongeLayer(const mu::Parser& parser);
    void setSpongeLayer(const std::string& muParserString);

    void setBlock(std::shared_ptr<Block3D> block);
    std::shared_ptr<Block3D> getBlock() const;

    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const;


    void swapDistributions();

protected:
    std::shared_ptr<DataSet3D> dataSet;
    std::shared_ptr<BCProcessor> bcProcessor;
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

    std::weak_ptr<Block3D> block;

    int nx1, nx2, nx3;

private:
    void checkFunction(mu::Parser fct);
};

#endif
