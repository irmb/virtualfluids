//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef SPHERE_MOCKS_H
#define SPHERE_MOCKS_H


#include "../Object.h"
#include <VirtualFluidsDefinitions.h>
#include <core/DataTypes.h>
#include "core/PointerDefinitions.h"

class SphereDummy : public Object
{
public:
    double getX1Centroid() override { return 0.0; }
    double getX1Minimum() override { return 0.0; }
    double getX1Maximum() override { return 0.0; }
    double getX2Centroid() override { return 0.0; }
    double getX2Minimum() override { return 0.0; }
    double getX2Maximum() override { return 0.0; }
    double getX3Centroid() override { return 0.0; }
    double getX3Minimum() override { return 0.0; }
    double getX3Maximum() override { return 0.0; }
    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset,
        const double& maxOffset) override {}
    bool isOnBoundary(const double& x1, const double& x2, const double& x3, const double& minOffset,
        const double& maxOffset) override {}
};

class SphereStub : public SphereDummy
{
private:
    SphereStub(double centerX, double centerY, double centerZ, double radius)
        : x(centerX), y(centerY), z(centerZ), radius(radius) {}
public:
    SPtr<SphereStub> makeShared(double centerX, double centerY, double centerZ, double radius) {
        return SPtr<SphereStub>(new SphereStub(centerX, centerY, centerZ, radius));
    }

    double getX1Minimum() override { return 0.0; }
    double getX1Maximum() override { return 0.0; }
    double getX2Minimum() override { return 0.0; }
    double getX2Maximum() override { return 0.0; }
    double getX3Minimum() override { return 0.0; }
    double getX3Maximum() override { return 0.0; }

private:
    double x, y, z;
    double radius;
};



#endif   
