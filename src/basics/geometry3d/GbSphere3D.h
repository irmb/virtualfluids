//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBSPHERE3D_H
#define GBSPHERE3D_H

#ifdef CAB_CTL
#include <ctl.h>
#endif // CAB_CTL

#include <cmath>
#include <vector>

#include <basics/utilities/UbObserver.h>

#include <geometry3d/GbObject3D.h>
#include <geometry3d/GbPoint3D.h>

#include <PointerDefinitions.h>

class GbLine3D;
class GbTriangle3D;
class GbObject3DCreator;

class GbSphere3D : public GbObject3D, public UbObserver
{
public:
    enum TRIANGULATIONMODE { CUBOIDPROJECTION, RAYPROJECTION };

    //////////////////////////////////////////////////////////////////////////
    // Konstruktoren
    GbSphere3D();
    GbSphere3D(const double &x1, const double &x2, const double &x3, const double &radius);
    GbSphere3D(const GbSphere3D &sphere);
    GbSphere3D(GbSphere3D *sphere); //<-unschoen!

    ~GbSphere3D() override;

    GbSphere3D *clone() override { return new GbSphere3D(*this); }
    void finalize() override;

    bool intersects(SPtr<GbSphere3D> sphere);

    double getRadius() const { return this->radius; }

    double getX1Centroid() override { return midPoint->getX1Coordinate(); }
    double getX1Minimum() override { return midPoint->getX1Coordinate() - radius; }
    double getX1Maximum() override { return midPoint->getX1Coordinate() + radius; }
    double getX2Centroid() override { return midPoint->getX2Coordinate(); }
    double getX2Minimum() override { return midPoint->getX2Coordinate() - radius; }
    double getX2Maximum() override { return midPoint->getX2Coordinate() + radius; }
    double getX3Centroid() override { return midPoint->getX3Coordinate(); }
    double getX3Minimum() override { return midPoint->getX3Coordinate() - radius; }
    double getX3Maximum() override { return midPoint->getX3Coordinate() + radius; }

    void setCenterX1Coordinate(const double &value) override;
    void setCenterX2Coordinate(const double &value) override;
    void setCenterX3Coordinate(const double &value) override;
    void setCenterCoordinates(const double &x1, const double &x2, const double &x3) override;
    void setCenterCoordinates(const UbTupleDouble3 &position) override;
    void setRadius(const double &radius);

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;
    double getDistance(GbPoint3D *p);
    double getDistance(const double &x1p, const double &x2p, const double &x3p);

    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3) override;
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary) override;

    bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                 const double &x2b, const double &x3b) override;
    bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;
    double getCellVolumeInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;
    double getCellVolumeInsideGbObject3DHelperFunction(const double &x1a, const double &x2a, const double &x3a,
                                                       const double &x1b, const double &x2b, const double &x3b);

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override;

    bool hasRaytracing() override { return true; }
    /*|r| must be 1! einheitsvector!!*/
    double getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3, const double &rx1,
                                         const double &rx2, const double &rx3) override;

    bool hasIntersectionWithDirectedLine(GbPoint3D origin, GbPoint3D direction);

    std::string toString() override;

    void translate(const double &x1, const double &x2, const double &x3) override
    {
        this->midPoint->translate(x1, x2, x3);
        this->notifyObserversObjectChanged();
    }
    void rotate(const double &rx1, const double &rx2, const double &rx3) override
    { /* rotation makes no sense*/
    }
    void scale(const double &sx1, const double & /*sx2*/, const double & /*sx3*/) override
    {
        this->radius *= sx1;
        this->notifyObserversObjectChanged();
    }

    void transform(const double matrix[4][4]);

    TRIANGULATIONMODE getTriangulationMode() { return triangulationMode; }
    void setTriangulationMode(TRIANGULATIONMODE mode) { this->triangulationMode = mode; }

    // virtuelle Methoden von UbObserver
    void objectChanged(UbObservable * /*changedObject*/) override
    {
        this->notifyObserversObjectChanged();
        // std::cout<<"GbSphere:objectChanged() - toDo-);";
    }
    void objectWillBeDeleted(UbObservable * /*objectForDeletion*/) override
    {
        throw UbException(UB_EXARGS, "not implemented");
    }

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere, weil man eine

#ifdef CAB_CTL
    ctl::oStream &write(ctl::oStream &os) const
    {
        midPoint->write(os);
        return os << radius;
    }
    ctl::iStream &read(ctl::iStream &is)
    {
        midPoint->read(is);
        return is >> radius;
    }
#endif // CAB_CTL

private:
    GbPoint3D *midPoint;
    double radius; // Radius des Kreises
    TRIANGULATIONMODE triangulationMode;
};

#endif // GBSPHERE3D_H
