//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef OBJECT_H
#define OBJECT_H

#include <VirtualFluidsDefinitions.h>
#include "grid/Cell.h"
#include "global.h"

class GridImp;
struct Vertex;

class VF_PUBLIC Object
{
public:
    HOSTDEVICE virtual ~Object() {}
    HOSTDEVICE virtual Object* clone() const = 0;

    virtual double getX1Centroid() = 0;
    virtual double getX1Minimum()  = 0;
    virtual double getX1Maximum()  = 0;

    virtual double getX2Centroid() = 0;
    virtual double getX2Minimum()  = 0;
    virtual double getX2Maximum()  = 0;

    virtual double getX3Centroid() = 0;
    virtual double getX3Minimum()  = 0;
    virtual double getX3Maximum()  = 0;


    virtual void scale(double delta) = 0;


    HOSTDEVICE virtual bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) = 0;

    HOSTDEVICE virtual bool isCellInObject(const Cell& cell) {
        for (const auto point : cell)
        {
            const bool isInObject = isPointInObject(point.x, point.y, point.z, 0.0, 0.0);
            if (!isInObject)
                return false;
        }
        return true;
    }

    CUDA_HOST virtual void findInnerNodes(SPtr<GridImp> grid);

    CUDA_HOST virtual int getIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnObject, real &qVal);
};


#endif
