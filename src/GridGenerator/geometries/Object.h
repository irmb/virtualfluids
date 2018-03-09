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


class VF_PUBLIC Object
{
public:
    HOSTDEVICE virtual ~Object() {};

    virtual double getX1Centroid() = 0;
     double getX1Minimum() { return x1min; }
     double getX1Maximum() { return x1max; }

    virtual double getX2Centroid() = 0;
     double getX2Minimum() { return x2min; }
     double getX2Maximum() { return x2max; }

    virtual double getX3Centroid() = 0;
     double getX3Minimum() { return x3min; }
     double getX3Maximum() { return x3max; }

    virtual void setX1Minimum(double value) { this->x1min = value; }
    virtual void setX1Maximum(double value) { this->x1max = value; }
    
    virtual void setX2Minimum(double value) { this->x2min = value; }
    virtual void setX2Maximum(double value) { this->x2max = value; }
    
    virtual void setX3Minimum(double value) { this->x3min = value; }
    virtual void setX3Maximum(double value) { this->x3max = value; }

    ///*=======================================================*/
    //double getLengthX1() { return (getX1Maximum() - getX1Minimum()); }
    //double getLengthX2() { return (getX2Maximum() - getX2Minimum()); }
    //double getLengthX3() { return (getX3Maximum() - getX3Minimum()); }

    //virtual void setCenterX1Coordinate(const double& value) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void setCenterX2Coordinate(const double& value) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void setCenterX3Coordinate(const double& value) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void setCenterCoordinates(const double& x1, const double& x2, const double& x3) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void setCenterCoordinates(const UbTupleDouble3& position) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }

    ////Rotates the Point in relation to the origin.
    ////Parameters must be radian measure.
    //virtual void rotate(const double& rx1, const double& rx2, const double& rx3) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void translate(const double& x1, const double& x2, const double& x3) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
    //virtual void scale(const double& sx1, const double& sx2, const double& sx3) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }


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

    virtual bool isOnBoundary(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) = 0;

    //virtual bool isPointInObject(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary) = 0;
    //virtual bool isPointInObject(const double& x1, const double& x2, const double& x3) = 0;

    //virtual bool isCellInsideGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b);
    //virtual bool isCellCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b);
    //virtual bool isCellInsideOrCuttingGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b);
    //virtual double getCellVolumeInsideGbObject3D(const double& x1a, const double& x2a, const double& x3a, const double& x1b, const double& x2b, const double& x3b) { return -1.0; };

    //virtual bool isInsideCell(const double& minX1, const double& minX2, const double& minX3, const double& maxX1, const double& maxX2, const double& maxX3);

protected:
    double x1min, x2min, x3min, x1max, x2max, x3max;

};


#endif
