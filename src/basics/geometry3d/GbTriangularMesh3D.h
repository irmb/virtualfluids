//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBTRIANGULARMESH_H
#define GBTRIANGULARMESH_H

#include <iostream>
#include <sstream>

#include <geometry3d/GbLine3D.h>
#include <geometry3d/GbPoint3D.h>
#include <geometry3d/GbTriangle3D.h>
#include <geometry3d/GbVector3D.h>

#include <basics/writer/WbWriter.h>

#ifdef CAB_RCF
#include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif // CAB_RCF

#include <PointerDefinitions.h>

/*=========================================================================*/
/* GbTriangularMesh3D                                                                  */
/*                                                                         */
/**
 * This Class provides the triangular meshes.
 * Note, that up to now no methods for checking consistency are included.
 * in this context this class describes facettes from an 3D-object !!!
 */
class GbTriangularMesh3D : public GbObject3D
{
public:
    enum POINTINOBJECTTEST { RAYCROSSING, HALFSPACE };

    GbTriangularMesh3D();
    GbTriangularMesh3D(std::string name, std::vector<GbPoint3D *> *nodes, std::vector<GbTriangle3D *> *triangles);
    GbTriangularMesh3D(std::string name, std::vector<GbTriangle3D *> *triangles);
    GbTriangularMesh3D(std::string name, std::vector<GbPoint3D *> *nodes, std::vector<GbLine3D *> *edges,
                       std::vector<GbTriangle3D *> *triangles);
    ~GbTriangularMesh3D() override;
    GbTriangularMesh3D *clone() override { throw UbException(UB_EXARGS, "not implemented"); }
    void finalize() override { throw UbException("GbTriangularMesh3D::finalize() - toDo"); }
    void setPointInObjectTest(POINTINOBJECTTEST mode) { this->pointinobjecttest = mode; }

    std::string getClassName() { return "GbTriangularMesh3D"; }

    std::string toString() override;
    // std::string getName();
    std::vector<GbPoint3D *> *getNodes();
    std::vector<GbTriangle3D *> *getTriangles();
    double getX1Centroid() override;
    double getX2Centroid() override;
    double getX3Centroid() override;
    double getX1Minimum() override;
    double getX1Maximum() override;
    double getX2Minimum() override;
    double getX2Maximum() override;
    double getX3Minimum() override;
    double getX3Maximum() override;

    void rotate(const double &alpha, const double &beta, const double &gamma) override;
    void translate(const double &x1, const double &x2, const double &x3) override;

    void calculateValues();
    void deleteRedundantNodes();

    UbTupleDouble6 calculateMomentOfInertia(double rhoP);
    UbTupleDouble3 calculateCenterOfGravity();

    double getArea();
    double getVolume();
    double getVolumeForRectangle(const double &p1x1, const double &p1x2, const double &p2x1, const double &p2x2);
    std::vector<GbTriangle3D *> *getTrianglesForRectangle(const double &p1x1, const double &p1x2, const double &p2x1,
                                                          const double &p2x2);
    std::vector<GbPoint3D *> *getNodesForRectangle(const double &p1x1, const double &p1x2, const double &p2x1,
                                                   const double &p2x2);
    double getX3RangeForRectangle(const double &p1x1, const double &p1x2, const double &p2x1, const double &p2x2);
    double getX3MinimumForRectangle(const double &p1x1, const double &p1x2, const double &p2x1, const double &p2x2);
    double getX3MaximumForRectangle(const double &p1x1, const double &p1x2, const double &p2x1, const double &p2x2);

    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3) override;

    bool isPointInGbObject3D(const double &x1, const double &x2, const double &x3, bool &pointIsOnBoundary) override;

    bool isPointInObject3DHalfSpace(const double &xp, const double &yp,
                                    const double &zp); // based on Halfspace algorithm
    bool isPointInObject3DSpherical(const double &xp, const double &yp, const double &zp,
                                    int numTriangles); // based on Spherical polygon area method

    // should be checked !!!
    bool isPointInObject3DRayCrossing(const double &xp, const double &yp, const double &zp, int radius, int numVertices,
                                      int numTriangles); // based on Ray tracing algorithm

    bool InPolyhedron(int F, GbVector3D &q, int radius);
    void RandomRay(GbVector3D &ray, int radius);
    char SegPlaneInt(GbTriangle3D *Tri, GbVector3D &q, GbVector3D &r, GbVector3D &p, int *m);
    int PlaneCoeff(GbTriangle3D *Tri, GbVector3D &Normal, double *D);
    char InTri3D(GbTriangle3D *T, int m, GbVector3D &p);
    char InTri2D(GbVector3D Tp[3], GbVector3D &pp);
    double AreaSign(GbVector3D &a, GbVector3D &b, GbVector3D &c);
    char SegTriInt(GbTriangle3D *Tri, GbVector3D &q, GbVector3D &r, GbVector3D &p);
    char InPlane(GbTriangle3D *T, int m, GbVector3D &q, GbVector3D &r, GbVector3D &p);
    char SegTriCross(GbTriangle3D *T, GbVector3D &q, GbVector3D &r);
    double VolumeSign(GbVector3D &a, GbVector3D &b, GbVector3D &c, GbVector3D &d);
    bool BoxTest(GbTriangle3D *triangle, GbVector3D &PointQ, GbVector3D &PointR);
    // till here !!!

    GbLine3D *createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2) override;
    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;

    void writeMesh(std::string filename, WbWriter *writer, bool writeNormals = false);
    /*======================================================================*/
    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

#ifdef CAB_RCF
    template <class Archive>
    void SF_SERIALIZE(Archive &ar)
    {
        SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
        ar &triangles;
        if (ArchiveTools::isWriting(ar)) {
            for (std::size_t t = 0; t < triangles->size(); t++) {
                nodes->push_back((*triangles)[t]->getPoint(0));
                nodes->push_back((*triangles)[t]->getPoint(1));
                nodes->push_back((*triangles)[t]->getPoint(2));
            }
        }
        // ar & nodes; //<- problem redundanz
        // ar & edges;
        ar &pointinobjecttest;
        ar &x1min;
        ar &x1max;
        ar &x2min;
        ar &x2max;
        ar &x3min;
        ar &x3max;
        ar &consistent;
    }
#endif // CAB_RCF

protected:
    std::vector<GbPoint3D *> *nodes;
    std::vector<GbLine3D *> *edges;
    std::vector<GbTriangle3D *> *triangles;

private:
    POINTINOBJECTTEST pointinobjecttest;
    void init();
    /*======================================================================*/
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double x3min;
    double x3max;
    bool consistent;
};
/*=========================================================================*/

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
#if CAB_RCF <= 903
SF_SERIALIZE_ENUM(GbTriangularMesh3D::POINTINOBJECTTEST) // bei klassen ausserhalb der klasse;-)
#endif
UB_AUTO_RUN_NAMED(SF::registerType<GbTriangularMesh3D>("GbTriangularMesh3D "), SF_GbTriangularMesh3D);
UB_AUTO_RUN_NAMED((SF::registerBaseAndDerived<GbObject3D, GbTriangularMesh3D>()), SF_GbTriangularMesh3D_BD1);
#endif // RCF_USE_SF_SERIALIZATION

#endif
