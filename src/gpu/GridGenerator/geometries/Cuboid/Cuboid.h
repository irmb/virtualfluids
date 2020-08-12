//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CUBOID_H
#define CUBOID_H

#include "global.h"
#include "GridGenerator_export.h"

#include "geometries/Object.h"

class GRIDGENERATOR_EXPORT Cuboid : public Object
{
public:              
    HOSTDEVICE Cuboid(const double& minX1, const double& minX2, const double& minX3, const double& maxX1,const double& maxX2, const double& maxX3);
    HOSTDEVICE virtual ~Cuboid();

    HOSTDEVICE Object* clone() const override;

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

    HOSTDEVICE bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;

private:
    static double getCenter(double x1, double x2);
    static double getMinimum(double x1, double x2);
    static double getMaximum(double x1, double x2);
    static bool isOn(const real& coord, const real& plane1, const real& plane2);
    static bool isBetween(const real& coord, const real& start, const real& end);

protected:
    double minX1;
    double minX2;
    double minX3;
    double maxX1;
    double maxX2;
    double maxX3;
};



#endif   
