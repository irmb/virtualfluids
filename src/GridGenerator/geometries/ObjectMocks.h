#ifndef OBJECT_MOCKS_H
#define OBJECT_MOCKS_H

#include "GridGenerator/global.h"

#include "geometries/Object.h"

class ObjectDummy : public Object
{
public:
    virtual ~ObjectDummy() {}

    static SPtr<ObjectDummy> makeShared()
    {
        return SPtr<ObjectDummy>(new ObjectDummy());
    }

    virtual Object* clone() const override { return nullptr; }
    virtual double getX1Centroid() override { return 0.0; }
    virtual double getX1Minimum() override { return 0.0; }
    virtual double getX1Maximum() override { return 0.0; }
    virtual double getX2Centroid() override { return 0.0; }
    virtual double getX2Minimum() override { return 0.0; }
    virtual double getX2Maximum() override { return 0.0; }
    virtual double getX3Centroid() override { return 0.0; }
    virtual double getX3Minimum() override { return 0.0; }
    virtual double getX3Maximum() override { return 0.0; }
    virtual void scale(double delta) override {}
    virtual bool isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset) override { return false; }
};


class ObjectStub : public ObjectDummy
{
public:
    virtual ~ObjectStub() {}

    static SPtr<ObjectStub> makeShared()
    {
        return SPtr<ObjectStub>(new ObjectStub());
    }
    double getX1Minimum() override {
        return this->minX1;
    }
    double getX1Maximum() override { return this->minX2; }
    double getX2Minimum() override { return this->minX3; }
    double getX2Maximum() override { return this->maxX1; }
    double getX3Minimum() override { return this->maxX2; }
    double getX3Maximum() override { return this->maxX3; }

    void setX1Minimum(double minX1) { this->minX1 = minX1; }
    void setX1Maximum(double minX2) { this->minX2 = minX2; }
    void setX2Minimum(double minX3) { this->minX3 = minX3; }
    void setX2Maximum(double maxX1) { this->maxX1 = maxX1; }
    void setX3Minimum(double maxX2) { this->maxX2 = maxX2; }
    void setX3Maximum(double maxX3) { this->maxX3 = maxX3; }


protected:
    double minX1, minX2, minX3, maxX1, maxX2, maxX3;
};




#endif
