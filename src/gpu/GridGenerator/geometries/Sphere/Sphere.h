//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef SPHERE_H
#define SPHERE_H

#include "global.h"

#include "geometries/Object.h"

class VF_PUBLIC Sphere : public Object
{
public:
    HOSTDEVICE Sphere(const double& centerX, const double& centerY, const double& centerZ, const double& radius);
    HOSTDEVICE virtual ~Sphere();

    static SPtr<Sphere> makeShared(double centerX, double centerY, double centerZ, double radius);

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
    
    CUDA_HOST int getIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnObject, real &qVal) override;


protected:
    double centerX;
    double centerY;
    double centerZ;

    double radius;
};



#endif   
