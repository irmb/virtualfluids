#ifndef KDTREE_H
#define KDTREE_H

#include <basics/utilities/UbKeys.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/KdNode.h>
#include <geometry3d/KdTree/KdRay.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>
#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

#include <cmath>
#include <string>

namespace Kd
{
template <typename T>
class Tree
{
public:
    /* ======================================================================================= */
    Tree(GbTriFaceMesh3D &mesh, const SplitAlgorithm<T> &splitAlg) : rootNode(NULL) { this->buildTree(mesh, splitAlg); }
    /* ======================================================================================= */
    ~Tree()
    {
        if (rootNode) {
            delete rootNode;
            rootNode = NULL;
        }
    }
    /* ======================================================================================= */
    // the IntersectionHandler specifies how to handle the intersection
    bool intersectLine(const UbTuple<T, T, T> &n1, const UbTuple<T, T, T> &n2,
                       const LineIntersectionHandler<T> &iHandler)
    {
        return rootNode->intersectLine(n1, n2, iHandler);
    }
    /* ======================================================================================= */
    // the IntersectionHandler specifies how to handle the intersection
    int intersectRay(const Ray<T> &ray, const RayIntersectionHandler<T> &iHandler)
    {
        std::set<UbKeys::Key3<int>> mailbox;
        return rootNode->intersectRay(ray, iHandler, mailbox);
    }
    /* ======================================================================================= */
    int getNumOfNodes()
    {
        if (rootNode)
            return rootNode->getNumOfNodes();
        return 0;
    }
    /* ======================================================================================= */
    int getNumOfTriFaces()
    {
        if (rootNode)
            return rootNode->getNumOfTriFaces();
        return 0;
    }
    /* ======================================================================================= */
    std::string toString()
    {
        return ""; // Tree:: num of nodes: " + rootNode.getNumOfNodes() + ", primitives:" +
                   // rootNode.getNumOfPrimitives() + ", root_primitives:" + getNumOfPrimitives() + ", max_level:" +
                   // max_level;
    }
    /* ======================================================================================= */
    void buildTree(GbTriFaceMesh3D &mesh, const SplitAlgorithm<T> &splitAlg)
    {
        if (rootNode)
            delete rootNode;

        // create a copy of triangles
        MbSmartPtr<std::vector<GbTriFaceMesh3D::TriFace>> triFaces(
            new std::vector<GbTriFaceMesh3D::TriFace>(*mesh.getTriangles()));

        const int maxLevel =
            static_cast<int>(lround(8.0 + 1.3 * std::log((double)triFaces->size()))); // TODO: remove magic numbers

        rootNode =
            new Node<T>(T(mesh.getX1Minimum()), T(mesh.getX2Minimum()), T(mesh.getX3Minimum()), T(mesh.getX1Maximum()),
                        T(mesh.getX2Maximum()), T(mesh.getX3Maximum()), triFaces, mesh.getNodes());

        rootNode->buildTree(0, maxLevel, splitAlg);
    }
    void writeTree(const std::string &filename, WbWriter *writer = WbWriterVtkXmlBinary::getInstance())
    {
        if (rootNode) {
            std::vector<UbTupleFloat3> nodes;
            std::vector<UbTupleInt8> cubes;
            std::vector<std::string> datanames;
            std::vector<std::vector<double>> cubesdata;
            rootNode->addCubeInfo(nodes, cubes, datanames, cubesdata);
            writer->writeOctsWithCellData(filename, nodes, cubes, datanames, cubesdata);
        }
    }

private:
    Node<T> *rootNode;
};
} // namespace Kd

#endif // KDTREE_H
