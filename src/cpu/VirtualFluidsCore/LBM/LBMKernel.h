#ifndef LBMKERNEL_H
#define LBMKERNEL_H

#include <PointerDefinitions.h>
#include "LBMSystem.h"
#include "ILBMKernel.h"
#include <array>
#include <muParser.h>

class BCProcessor;
class DataSet3D;
class Block3D;

class LBMKernel : public ILBMKernel, public enableSharedFromThis<LBMKernel>
{
public:
    typedef std::numeric_limits<LBMReal> LBMRealLim;
public:
    LBMKernel();
    virtual ~LBMKernel();

    virtual SPtr<LBMKernel> clone() = 0;

    virtual void calculate(int step) = 0;
    virtual double getCalculationTime() = 0;

    void setBCProcessor(SPtr<BCProcessor> bcp);
    SPtr<BCProcessor> getBCProcessor() const;

    void setCollisionFactor(double collFactor);
    double getCollisionFactor() const;

    void setGhostLayerWidth(int witdh);
    int  getGhostLayerWidth() const;

    void setDataSet(SPtr<DataSet3D> dataSet);
    SPtr<DataSet3D> getDataSet() const;

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

    void setBlock(SPtr<Block3D> block);
    SPtr<Block3D> getBlock() const;

    bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const;

    void swapDistributions();

    void setNX(std::array<int, 3> nx);
    std::array<int, 3> getNX();

protected:
    SPtr<DataSet3D> dataSet;
    SPtr<BCProcessor> bcProcessor;
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

    WPtr<Block3D> block;

    std::array<int, 3> nx;

private:
    void checkFunction(mu::Parser fct);
};

#endif
