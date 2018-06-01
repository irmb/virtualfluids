//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef VERTICAL_CYLINDER_H
#define VERTICAL_CYLINDER_H


#include "../Object.h"
#include <VirtualFluidsDefinitions.h>
#include <core/DataTypes.h>
#include <core/PointerDefinitions.h>

class VF_PUBLIC VerticalCylinder : public Object
{
public:
    HOSTDEVICE VerticalCylinder(const double& centerX, const double& centerY, const double& centerZ, const double& radius, const double& height);
    HOSTDEVICE virtual ~VerticalCylinder();

    static SPtr<VerticalCylinder> makeShared(double centerX, double centerY, double centerZ, double radius, double height);

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

    HOSTDEVICE bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;


    void scale(double delta) override;
   
protected:
    double centerX;
    double centerY;
    double centerZ;

    double radius;
    double height;
};



#endif   
