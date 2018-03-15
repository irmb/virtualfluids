//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef CONGLOMERATE_H
#define CONGLOMERATE_H


#include "../Object.h"
#include <VirtualFluidsDefinitions.h>
#include <core/DataTypes.h>
#include "global.h"

#define MAX_NUMBER_OF_OBJECTS 20

class VF_PUBLIC Conglomerate : public Object
{
public:              
    HOSTDEVICE Conglomerate();
    HOSTDEVICE virtual ~Conglomerate();

    HOST void add(Object* object);

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

protected:
    Object** objects;
    uint numberOfObjects = 0;
};



#endif   
