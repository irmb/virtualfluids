//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef SPHERE_H
#define SPHERE_H


#include "../Object.h"
#include <VirtualFluidsDefinitions.h>
#include <core/DataTypes.h>
#include <core/PointerDefinitions.h>

class VF_PUBLIC Sphere : public Object
{
public:
    HOSTDEVICE Sphere(const double& centerX, const double& centerY, const double& centerZ, const double& radius);
    HOSTDEVICE virtual ~Sphere();

    static SPtr<Sphere> makeShared(double centerX, double centerY, double centerZ, double radius);

    HOSTDEVICE Object* clone() const;

    double getX1Centroid() override;
    double getX1Minimum();
    double getX1Maximum();
    double getX2Centroid() override;
    double getX2Minimum();
    double getX2Maximum();
    double getX3Centroid() override;
    double getX3Minimum();
    double getX3Maximum();

    HOSTDEVICE bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;
    bool isOnBoundary(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;


public:
/*
    void translate(const double& x1, const double& x2, const double& x3);
    void scale(const double& sx1, const double& sx2, const double& sx3);

    double getLengthX1();
    double getLengthX2();
    double getLengthX3();

    bool isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
    bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
    bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
    double getCellVolumeInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);


    bool isPointInObject(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary) override;
    bool isPointInObject(const double& x1, const double& x2, const double& x3) override;*/
protected:
    double centerX;
    double centerY;
    double centerZ;

    double radius;
};



#endif   
