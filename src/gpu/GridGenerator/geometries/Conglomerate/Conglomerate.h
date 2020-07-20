//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CONGLOMERATE_H
#define CONGLOMERATE_H

#include "global.h"

#include "geometries/Object.h"

#define MAX_NUMBER_OF_OBJECTS 20

class VIRTUALFLUIDS_GPU_EXPORT Conglomerate : public Object
{
public:              
    HOSTDEVICE Conglomerate();
    HOSTDEVICE virtual ~Conglomerate();

    CUDA_HOST static SPtr<Conglomerate> makeShared();

    HOSTDEVICE void add(Object* object);
    HOSTDEVICE void subtract(Object* objectStub);


    HOSTDEVICE    Object* clone() const override;

    double getX1Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Centroid() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Centroid() override;
    double getX3Minimum() override;
    double getX3Maximum() override;

    void scale(double delta) override;

    HOSTDEVICE    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;

    CUDA_HOST void findInnerNodes(SPtr<GridImp> grid) override;

protected:
    static double getMinimum(double val1, double val2);
    static double getMaximum(double val1, double val2);


    Object** addObjects;
    Object** subtractObjects;
    uint numberOfAddObjects = 0;
    uint numberOfSubtractObjects = 0;
};



#endif   
