//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CUBOID_H
#define CUBOID_H


#include "../Object.h"
#include <VirtualFluidsDefinitions.h>
#include <core/DataTypes.h>

class VF_PUBLIC Cuboid : public Object
{
public:              
    Cuboid(const double& minX1, const double& minX2, const double& minX3, const double& maxX1,const double& maxX2, const double& maxX3);
    virtual ~Cuboid();

    Object* clone() const;

    double getX1Centroid() override;
    double getX1Minimum();
    double getX1Maximum();
    double getX2Centroid() override;
    double getX2Minimum();
    double getX2Maximum();
    double getX3Centroid() override;
    double getX3Minimum();
    double getX3Maximum();

    bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;

    bool isOnBoundary(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override;

private:
    static double getCenter(double x1, double x2);
    static double getMinimum(double x1, double x2);
    static double getMaximum(double x1, double x2);
    static bool isOn(const real& coord, const real& plane1, const real& plane2);
    static bool isBetween(const real& coord, const real& start, const real& end);

public:
    /*
       double getX1Centroid();
       double getX1Minimum();
       double getX1Maximum();
       double getX2Centroid();
       double getX2Minimum();
       double getX2Maximum();
       double getX3Centroid();
       double getX3Minimum();
       double getX3Maximum();

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
    double minX1;
    double minX2;
    double minX3;
    double maxX1;
    double maxX2;
    double maxX3;
};



#endif   
