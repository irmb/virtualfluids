#ifndef I_LBMKERNEL_H
#define I_LBMKERNEL_H

#include <PointerDefinitions.h>


class BCProcessor;
class DataSet3D;

class ILBMKernel
{
public:
    virtual ~ILBMKernel() {};

    virtual void calculate(int step) = 0;
    virtual double getCalculationTime() = 0;
    virtual void swapDistributions() = 0;

    virtual bool getCompressible() const = 0;
    virtual SPtr<BCProcessor> getBCProcessor() const = 0;
    virtual void setBCProcessor(SPtr<BCProcessor> bcProcessor) = 0;
    virtual SPtr<DataSet3D> getDataSet() const = 0;
    virtual double getCollisionFactor() const = 0;
    virtual void setCollisionFactor(double collFactor) = 0;
    virtual bool isInsideOfDomain(const int &x1, const int &x2, const int &x3) const = 0;
    virtual int getGhostLayerWidth() const = 0;
    virtual double getDeltaT() const = 0;
    virtual bool getWithForcing() const = 0;
};

#endif
