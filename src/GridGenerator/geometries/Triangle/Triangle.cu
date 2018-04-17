#include "Triangle.h"
#include <GridGenerator/utilities/math/Math.h>
#include "TriangleException.h"

#include "Serialization/TriangleMemento.h"

#include "grid/NodeValues.h"

HOSTDEVICE Triangle::Triangle(Vertex &v1, Vertex &v2, Vertex &v3, Vertex &normal) : v1(v1), v2(v2), v3(v3), normal(normal) {}
HOSTDEVICE Triangle::Triangle(Vertex &v1, Vertex &v2, Vertex &v3) : v1(v1), v2(v2), v3(v3) { calcNormal(); }
HOSTDEVICE Triangle::Triangle(){}

HOSTDEVICE void Triangle::set(const Vertex &v1, const Vertex &v2, const Vertex &v3)
{
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->calcNormal();
}

HOSTDEVICE void Triangle::calcNormal()
{
    Vertex edge1 = v2 - v1;
    Vertex edge2 = v3 - v1;
    normal = edge1.crossProduct(edge2);
    normal.normalize();
}

HOSTDEVICE void Triangle::initalLayerThickness(real delta)
{
    this->layerThickness = delta*(abs(this->normal.x) + abs(this->normal.y) + abs(this->normal.z));
}


HOSTDEVICE char Triangle::isUnderFace(const Vertex &point) const
{
    real s;

    if (this->isUnterExtendedFace(point, s))
        if (this->isNotNextToFace(point))
            if (this->isUnderAngleToNeighbors(point))
                if (this->isNegativeDirectionBorder(point))
                    return NEGATIVE_DIRECTION_BORDER;
                else
                    return INSIDE;
            else 
                return 4;
        else
            return 3;


    if (this->isQNode(point, s))
        return Q;
    
    return FLUID;
}

HOSTDEVICE bool Triangle::isUnterExtendedFace(const Vertex & point, real &s) const
{
    s = this->getPerpedicularDistanceFrom(point);
    return ((vf::Math::greaterEqual(s, 0.0f)) && s < this->layerThickness);
}

HOSTDEVICE real Triangle::getPerpedicularDistanceFrom(const Vertex &P) const
{
    Vertex v = P - v1;
    return  (v * -1.0f) * normal;
}

HOSTDEVICE Vertex Triangle::getPerpedicularPointFrom(const Vertex &P) const
{
    return P + normal * getPerpedicularDistanceFrom(P);
}

HOSTDEVICE bool Triangle::isQNode(const Vertex & point, const real &s) const
{
    return (s < 0 && vf::Math::lessEqual(-s, this->layerThickness));
    //calculateQs(actualPoint, actualTriangle);
}

HOSTDEVICE bool Triangle::isNegativeDirectionBorder(const Vertex &point) const
{
    return normal.x < 0.0f || normal.y < 0.0f || normal.z < 0.0f;
    //return (sVector.x < 0.0f && sVector.y < 0.0f && sVector.z < 0.0f);
}

HOSTDEVICE bool Triangle::isNotNextToFace(const Vertex &point) const
{
    Vertex Pb = getPerpedicularPointFrom(point);

    Vertex w1 = Pb - v1;
    Vertex w2 = Pb - v2;
    Vertex w3 = Pb - v3;

    Vertex t1 = w1.crossProduct(v2 - v1);
    Vertex t2 = w2.crossProduct(v3 - v2);
    Vertex t3 = w3.crossProduct(v1 - v3);

    real g1 = t1 * normal;
    real g2 = t2 * normal;
    real g3 = t3 * normal;

    return vf::Math::lessEqual(g1, 0.0f) && vf::Math::lessEqual(g2, 0.0f) && vf::Math::lessEqual(g3, 0.0f);
}

HOSTDEVICE bool Triangle::isUnderAngleToNeighbors(const Vertex &point) const
{
    Vertex Pci[3];
    this->getClosestPointsOnEdges(Pci, point);
    Vertex Pb = this->getPerpedicularPointFrom(point);

    Vertex q[3];
    Vertex r[3];

    real betaAngles[3];
    for (int i = 0; i < 3; i++)
    {
        q[i] = point - Pci[i];
        r[i] = Pb - Pci[i];
        betaAngles[i] = q[i].getInnerAngle(r[i]);
    }

    real eps = EPSILON * 100.0f;
    return (vf::Math::lessEqual(betaAngles[0], alphaAngles[0], eps) && vf::Math::lessEqual(betaAngles[1], alphaAngles[1], eps) && vf::Math::lessEqual(betaAngles[2], alphaAngles[2], eps));
}

HOSTDEVICE void Triangle::getClosestPointsOnEdges(Vertex arr[], const Vertex &P) const
{
    Vertex Pc1, Pc2, Pc3;
    Vertex v4 = P - v1;
    Vertex v5 = P - v2;
    Vertex v6 = P - v3;

    Vertex d1 = v2 - v1;
    Vertex d2 = v3 - v2;
    Vertex d3 = v1 - v3;

    real temp = (v4 * d1) / (d1 * d1);
    Vertex tempV = d1 * temp;
    Pc1 = v1 + tempV;

    temp = (v5 * d2) / (d2 * d2);
    tempV = d2 * temp;
    Pc2 = v2 + tempV;

    temp = (v6 * d3) / (d3 * d3);
    tempV = d3 * temp;
    Pc3 = v3 + tempV;

    arr[0] = Pc1;
    arr[1] = Pc2;
    arr[2] = Pc3;
}

Vertex Triangle::getCenterOfMass() const 
{
    return (v1 + v2 + v3) * (1.0f / 3.0f);
}

real Triangle::getHalfAngleBetweenToAdjacentTriangle(const Triangle &t2) const
{
    if (isEqual(t2)) return 0.0f;

    real alpha = normal.getInnerAngle(t2.normal);
	if (alpha == 0.0f)
		return 90.0f;

    if(doesNormalsShowToEachOther(t2))
        return (180.0f + alpha) / 2.0f;
    else
        return (180.0f - alpha) / 2.0f;
}

int Triangle::isEqual(const Triangle &t2) const
{
    return getNumberOfCommonEdge(t2) == 3;
}

bool Triangle::doesNormalsShowToEachOther(const  Triangle &t2) const
{
    Vertex s1 = getCenterOfMass();
    Vertex s2 = t2.getCenterOfMass();

    Vertex s1s2 = s1 - s2;
    real X = s1s2 * t2.normal;
    return X > 0 ? true : false;
}

int Triangle::getCommonEdge(const Triangle &t2) const 
{
    bool edgeOneCommon = false;
	bool edgeTwoCommon = false;
	bool edgeThreeCommon = false;

	edgeOneCommon = t2.contains(v1);
	edgeTwoCommon = t2.contains(v2);
	edgeThreeCommon = t2.contains(v3);

    if (edgeOneCommon && edgeTwoCommon)
        return 0;
    else if (edgeTwoCommon && edgeThreeCommon)
        return 1;
    else if (edgeThreeCommon && edgeOneCommon)
        return 2;
    else
        return -1;
}

bool Triangle::contains(const Vertex& v) const 
{
	return (v == v1 || v == v2 || v == v3);
}


HOSTDEVICE int Triangle::getNumberOfCommonEdge(const Triangle &t2) const
{
	int commonEdge = 0;
	if (t2.contains(v1))
		commonEdge++;
	if (t2.contains(v2))
		commonEdge++;
	if (t2.contains(v3))
		commonEdge++;

    if (commonEdge == 2 || commonEdge == 3) return commonEdge;
    return 0;
}


HOSTDEVICE int Triangle::getTriangleIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnTri, real &qVal) const
{
    Vertex edge1, edge2, tvec, pvec, qvec, tuv;
    float det, inv_det;

    edge1 = v2 - v1;
    edge2 = v3 - v1;

    pvec = direction.crossProduct(edge2);
    det = edge1 * pvec;

    if (det < EPSILON)
        return 0;

    inv_det = 1 / det;

    tvec = P - v1;
    tuv.y = (tvec * pvec) * inv_det;

    if (tuv.y < 0.0 || tuv.y > 1.0)
        return 0;

    qvec = tvec.crossProduct(edge1);
    tuv.z = (direction * qvec) * inv_det;

    if (tuv.z < 0.0 || (tuv.y + tuv.z) > 1.0)
        return 0;

    tuv.x = (edge2 * qvec) * inv_det;

    pointOnTri.x = (1 - tuv.y - tuv.z) * v1.x + tuv.y * v2.x + tuv.z * v3.x;
    pointOnTri.y = (1 - tuv.y - tuv.z) * v1.y + tuv.y * v2.y + tuv.z * v3.y;
    pointOnTri.z = (1 - tuv.y - tuv.z) * v1.z + tuv.y * v2.z + tuv.z * v3.z;

    qVal = tuv.x;

    return 1;
}

HOSTDEVICE void Triangle::print() const
{
    printf("v1: ");
    v1.print();
    printf("v2: ");
    v2.print();
    printf("v3: ");
    v3.print();
    printf("normal: ");
    normal.print();
}

HOST bool Triangle::operator==(const Triangle &t) const
{
    return v1 == t.v1 && v2 == t.v2 && v3 == t.v3
        && vf::Math::equal(alphaAngles[0], t.alphaAngles[0]) && vf::Math::equal(alphaAngles[1], t.alphaAngles[1]) && vf::Math::equal(alphaAngles[2], t.alphaAngles[2]);
}


HOSTDEVICE void Triangle::setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ) const
{
    Vertex::setMinMax(minX, maxX, minY, maxY, minZ, maxZ, v1, v2, v3);
}

HOST TriangleMemento Triangle::getState() const
{
    TriangleMemento memento;
    memento.v1 = v1.getState();
    memento.v2 = v2.getState();
    memento.v3 = v3.getState();
    memento.normal = normal.getState();
    
    memento.alphaAngles[0] = alphaAngles[0];
    memento.alphaAngles[1] = alphaAngles[1];
    memento.alphaAngles[2] = alphaAngles[2];

    return memento;
}

HOST void Triangle::setState(const TriangleMemento &memento)
{
    this->v1.setState(memento.v1);
    this->v2.setState(memento.v2);
    this->v3.setState(memento.v3);
    this->normal.setState(memento.normal);

    this->alphaAngles[0] = memento.alphaAngles[0];
    this->alphaAngles[1] = memento.alphaAngles[1];
    this->alphaAngles[2] = memento.alphaAngles[2];
}

