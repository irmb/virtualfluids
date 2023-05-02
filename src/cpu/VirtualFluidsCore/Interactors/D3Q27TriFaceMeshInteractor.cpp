#include "D3Q27TriFaceMeshInteractor.h"
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbMath.h>

#include "basics/writer/WbWriterVtkXmlASCII.h"
#include <basics/writer/WbWriterVtkASCII.h>
#include <basics/writer/WbWriterVtkBinary.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "VelocityBC.h"
#include "basics/utilities/UbTiming.h"
#include <geometry3d/GbCuboid3D.h>
#include <geometry3d/GbHalfSpace3D.h>
#include <geometry3d/GbMeshTools3D.h>
#include <geometry3d/GbSystem3D.h>

#include <geometry3d/GbTriFaceMesh3D.h>

//#include <omp.h>

#include "CoordinateTransformation3D.h"
#include <stack>

using namespace std;

D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor() : D3Q27Interactor() { this->stressMode = STRESSNORMAL; }
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<Grid3D> /*grid*/, std::string /*name*/)
{
    this->stressMode = STRESSNORMAL;
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid,
                                                       SPtr<BC> BC, int type)
    : D3Q27Interactor(triFaceMesh, grid, BC, type)
{
    this->stressMode = STRESSNORMAL;
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::D3Q27TriFaceMeshInteractor(SPtr<GbTriFaceMesh3D> triFaceMesh, SPtr<Grid3D> grid,
                                                       SPtr<BC> BC, int type, Interactor3D::Accuracy a)
    : D3Q27Interactor(triFaceMesh, grid, BC, type, a)
{
    this->stressMode = STRESSNORMAL;
}
//////////////////////////////////////////////////////////////////////////
D3Q27TriFaceMeshInteractor::~D3Q27TriFaceMeshInteractor() = default;
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::initInteractor(const real &timeStep)
{
    updateBlocks(); 
    setQs(timeStep);
}
//////////////////////////////////////////////////////////////////////////
bool D3Q27TriFaceMeshInteractor::setDifferencesToGbObject3D(const SPtr<Block3D> block)
{
    if (!block)
        return false;

    // UBLOG(logINFO, "D3Q27TriFaceMeshInteractor::setDifferencesToGbObject3D()");

    bcNodeIndicesMap[block] = set<std::vector<int>>();
    //   set< std::vector<int> >& transNodeIndices = bcNodeIndicesMap[block];
    solidNodeIndicesMap[block]         = set<UbTupleInt3>();
    set<UbTupleInt3> &solidNodeIndices = solidNodeIndicesMap[block];

    bool oneEntryGotBC = false; // ob ueberhaupt ein eintrag ein BC zugewiesen wurde
    //   bool gotQs         = false; //true, wenn "difference" gesetzt wurde
    SPtr<BoundaryConditions> bc;

    SPtr<ILBMKernel> kernel = block->getKernel();
    SPtr<BCArray3D> bcArray = kernel->getBCSet()->getBCArray();

    real internX1, internX2, internX3;

    int startIX1 = 0, startIX2 = 0, startIX3 = 0;
    int stopIX1 = (int)bcArray->getNX1(), stopIX2 = (int)bcArray->getNX2(), stopIX3 = (int)bcArray->getNX3();

    //   double         dx       = grid.lock()->getDeltaX(block);
    //   UbTupleDouble3 orgDelta = grid.lock()->getNodeOffset(block);

    bool pointOnBoundary = false;

    for (int ix3 = startIX3; ix3 < stopIX3; ix3++) {
        for (int ix2 = startIX2; ix2 < stopIX2; ix2++) {
            for (int ix1 = startIX1; ix1 < stopIX1; ix1++) {
                Vector3D coords = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                internX1        = coords[0];
                internX2        = coords[1];
                internX3        = coords[2];

                if (this->isSolid()) {
                    if (this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3)) {
                        if (bcArray->isFluid(ix1, ix2, ix3)) {
                            solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                            bcArray->setSolid(ix1, ix2, ix3);
                        }
                    }
                } else if (this->isInverseSolid()) {
                    // bei inverse solid sind alle Knoten AUSSERHALB und auf der boundary SOLID
                    if (!this->geoObject3D->isPointInGbObject3D(internX1, internX2, internX3, pointOnBoundary) ||
                        pointOnBoundary == true) {
                        if (bcArray->isFluid(ix1, ix2, ix3)) {
                            solidNodeIndices.insert(UbTupleInt3(ix1, ix2, ix3));
                            bcArray->setSolid(ix1, ix2, ix3);
                        }
                    }
                }
            }
        }
    }

    return oneEntryGotBC;
}
//////////////////////////////////////////////////////////////////////////
// E.F. /4/16/2013
void D3Q27TriFaceMeshInteractor::setQs(const real &timeStep)
{
    using namespace vf::lbm::dir;

    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - setQs start ");
    if (!this->grid.lock())
        throw UbException(UB_EXARGS, "ups, no grid.lock()!!");

    if (this->reinitWithStoredQsFlag && !bcNodeIndicesAndQsMap.empty()) {
        this->reinitWithStoredQs(timeStep);
        return;
    }

    GbTriFaceMesh3D *mesh = dynamic_cast<GbTriFaceMesh3D *>(this->geoObject3D.get());

    //////////////////////////////////////////////////////////////////////////
    // init bcs
    //////////////////////////////////////////////////////////////////////////
    int nofAdapter = (int)this->BCs.size();
    if (nofAdapter == 0)
        std::cout
            << "WARNING - D3Q27TriFaceMeshInteractor::initInteractor Warning - no nodeAdapter available for " /*<<this->getName()*/
            << std::endl;
    bool needTimeDependence = false;
    for (int pos = 0; pos < nofAdapter; ++pos) {
        this->BCs[pos]->init(this, timeStep);
        if (this->BCs[pos]->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();

    //////////////////////////////////////////////////////////////////////////
    // grid.lock() info
    //////////////////////////////////////////////////////////////////////////
    int coarsestInitLevel = grid.lock()->getCoarsestInitializedLevel();
    int finestInitLevel   = grid.lock()->getFinestInitializedLevel();

    UbTupleInt3 blocknx = grid.lock()->getBlockNX();
    int blocknx1        = val<1>(blocknx); // gilt fuer alle Level
    int blocknx2        = val<2>(blocknx); // gilt fuer alle Level
    int blocknx3        = val<3>(blocknx); // gilt fuer alle Level

    // grobe Blocklaengen
    SPtr<CoordinateTransformation3D> trafo = grid.lock()->getCoordinateTransformator();
    double cblockDeltaX1, cblockDeltaX2, cblockDeltaX3, delta;
    cblockDeltaX1 = cblockDeltaX2 = cblockDeltaX3 = delta = 1.0 / (double)(1 << coarsestInitLevel);
    if (trafo) {
        cblockDeltaX1 = trafo->getX1CoordinateScaling() * delta;
        cblockDeltaX2 = trafo->getX2CoordinateScaling() * delta;
        cblockDeltaX3 = trafo->getX3CoordinateScaling() * delta;
    }
    // levelspezifische blocklaengen und knotenabstaende
    std::vector<std::vector<double>> nodeDeltaToNeigh(finestInitLevel + 1);
    std::vector<float> deltaMinX1(finestInitLevel + 1), deltaMinX2(finestInitLevel + 1),
        deltaMinX3(finestInitLevel + 1);
    std::vector<float> deltaMaxX1(finestInitLevel + 1), deltaMaxX2(finestInitLevel + 1),
        deltaMaxX3(finestInitLevel + 1);

    // Im Boltzmankontext muss dx1==dx2==dx3 sein!!
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX2 / (double)blocknx2));
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX3 / (double)blocknx3));

    for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
        double nodeDeltaX1 = cblockDeltaX1 / (double)(blocknx1 * (1 << (level - coarsestInitLevel)));
        double nodeDeltaX2 = cblockDeltaX2 / (double)(blocknx2 * (1 << (level - coarsestInitLevel)));
        double nodeDeltaX3 = cblockDeltaX3 / (double)(blocknx3 * (1 << (level - coarsestInitLevel)));

        std::vector<double> distNeigh(D3Q27System::FENDDIR + 1, 0.0);
        D3Q27System::calcDistanceToNeighbors(distNeigh, nodeDeltaX1, nodeDeltaX2, nodeDeltaX3);
        // D3Q27System::calcDistanceToNeighbors(distNeigh, nodeDeltaX1);

        nodeDeltaToNeigh[level].resize(D3Q27System::ENDDIR + 1, 0.0);
        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            nodeDeltaToNeigh[level][fdir] = distNeigh[fdir];
        }

        // im gegensatz zum allg. Cell3DInteractor kann man hier auf max(0.02*blockDeltaX1[level],fabs(...)) verzichten
        // da dies nur für blockDeltaCalculator->getMinX1Delta(level)==0.0 benötigt wird. ist im D3Q19... aber nie so
        // Geller: kann man nicht diesen befuckten DeltaCalculator weglassen und hier einfach die Formel zum Delta
        // rechnen reinpacken SirAnn: klar, mann kann auch weißwuerste am Alex verkaufen... aber zum einen ist das Ding
        // dazu da
        // und zum anderen sollt eman mal überlegen: "Formel zum Delta rechnen"->man muss rechnen
        // blockDeltaCalculator->getMinX1Delta(level) -> ein geinlinter wert wird geholt

        // TODO: set 5.0 as variable parameter in constructor, default 2.0
        deltaMinX1[level] = (float)(5.0 * nodeDeltaX1); // kein minus da unten -deltaMin
        deltaMinX2[level] = (float)(5.0 * nodeDeltaX2);
        deltaMinX3[level] = (float)(5.0 * nodeDeltaX3);
        deltaMaxX1[level] = (float)(5.0 * nodeDeltaX1);
        deltaMaxX2[level] = (float)(5.0 * nodeDeltaX2);
        deltaMaxX3[level] = (float)(5.0 * nodeDeltaX3);
    }

    //////////////////////////////////////////////////////////////////////////
    // bounding cubes des TriFaceMesh ermitteln (pro level)
    //////////////////////////////////////////////////////////////////////////
    // min/max Werte des Dreiecksnetzes holen
    //   double geoMinX1(0.0), geoMinX2(0.0), geoMinX3(0.0), geoMaxX1(0.0), geoMaxX2(0.0), geoMaxX3(0.0);

    //   geoMinX1 = this->geoObject3D->getX1Minimum();  geoMaxX1 = this->geoObject3D->getX1Maximum();
    //   geoMinX2 = this->geoObject3D->getX2Minimum();  geoMaxX2 = this->geoObject3D->getX2Maximum();
    //   geoMinX3 = this->geoObject3D->getX3Minimum();  geoMaxX3 = this->geoObject3D->getX3Maximum();

    //////////////////////////////////////////////////////////////////////////
    // DREIECKE: q-Bestimmung
    //////////////////////////////////////////////////////////////////////////

    // notwendige variablen initialisieren (u.a. blockDeltas des groben levels)
    float triPoints[3][3];
    float vx1 = 0.0, vx2 = 0.0, vx3 = 0.0;
    unsigned counterTriBoxOverlap = 0, counterAABBTriFace = 0, counterHalfspace = 0, counterBilligOBB = 0;
    std::vector<GbTriFaceMesh3D::TriFace> &triangles = *mesh->getTriangles();
    std::vector<GbTriFaceMesh3D::Vertex> &nodes      = *mesh->getNodes();
    std::map<SPtr<Block3D>, std::set<UbTupleInt3>> tmpSolidNodesFromOtherInteractors;

    int onePercent = UbMath::integerRounding(triangles.size() * 0.01);
    if (onePercent == 0)
        onePercent = 1;
    UbTimer setQTimer;
    setQTimer.start();
    UBLOG(logDEBUG3, " - setQs for " << (int)triangles.size() << " triangles");

    //   bool solidFromOtherInteractor = false;
    float blockMinX[3], blockMaxX[3], boxCenter[3], halfBoxSize[3];

    for (size_t t = 0; t < triangles.size(); t++) {
        //////////////////////////////////////////////////////////////////////////
        // Halfspace zum Dreieck generieren und min/max des Dreiecks ermitteln
        //////////////////////////////////////////////////////////////////////////
        GbTriFaceMesh3D::TriFace &triangle = triangles[t];

        GbTriFaceMesh3D::Vertex &v1 = nodes[triangle.v1];
        GbTriFaceMesh3D::Vertex &v2 = nodes[triangle.v2];
        GbTriFaceMesh3D::Vertex &v3 = nodes[triangle.v3];

        if (this->isInverseSolid()) {
            triangle.nx *= (-1);
            triangle.ny *= (-1);
            triangle.nz *= (-1);
        }
        GbHalfSpace3D halfSpace(v1.x, v1.y, v1.z, triangle.nx, triangle.ny, triangle.nz);

        //////////////////////////////////////////////////////////////////////////
        // fuer GbMeshTools3D::triBoxOverlap
        //////////////////////////////////////////////////////////////////////////
        triPoints[0][0] = v1.x;
        triPoints[0][1] = v1.y;
        triPoints[0][2] = v1.z;
        triPoints[1][0] = v2.x;
        triPoints[1][1] = v2.y;
        triPoints[1][2] = v2.z;
        triPoints[2][0] = v3.x;
        triPoints[2][1] = v3.y;
        triPoints[2][2] = v3.z;

        double minX1 = triangle.getMinX(nodes);
        double maxX1 = triangle.getMaxX(nodes);
        double minX2 = triangle.getMinY(nodes);
        double maxX2 = triangle.getMaxY(nodes);
        double minX3 = triangle.getMinZ(nodes);
        double maxX3 = triangle.getMaxZ(nodes);

        //////////////////////////////////////////////////////////////////////////
        // Schleife ueber alle Level
        //////////////////////////////////////////////////////////////////////////
        double e1x1, e1x2, e1x3, e2x1, e2x2, e2x3, px1, px2, px3, a, f, sx1, sx2, sx3, u, qx1, qx2, qx3, v;
        bool gotQs = false;
        SPtr<BoundaryConditions> bc;

        for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
            //////////////////////////////////////////////////////////////////////////
            // levelspezifisches BoundCube des Dreicks ermitteln und zugehörige Bloecke beziehen
            //////////////////////////////////////////////////////////////////////////
            double boundCubeTriangleMinX1 = minX1 - deltaMinX1[level];
            double boundCubeTriangleMaxX1 = maxX1 + deltaMaxX1[level];
            double boundCubeTriangleMinX2 = minX2 - deltaMinX2[level];
            double boundCubeTriangleMaxX2 = maxX2 + deltaMaxX2[level];
            double boundCubeTriangleMinX3 = minX3 - deltaMinX3[level];
            double boundCubeTriangleMaxX3 = maxX3 + deltaMaxX3[level];

            GbCuboid3D boundingCubeTriangle(boundCubeTriangleMinX1, boundCubeTriangleMinX2, boundCubeTriangleMinX3,
                                            boundCubeTriangleMaxX1, boundCubeTriangleMaxX2, boundCubeTriangleMaxX3);

            std::vector<SPtr<Block3D>> triBlocks;
            grid.lock()->getBlocksByCuboid(level, boundCubeTriangleMinX1, boundCubeTriangleMinX2,
                                           boundCubeTriangleMinX3, boundCubeTriangleMaxX1, boundCubeTriangleMaxX2,
                                           boundCubeTriangleMaxX3, triBlocks);

            //////////////////////////////////////////////////////////////////////////
            // Schleife ueber bloecke des level, die das dreieck beinhalten
            //////////////////////////////////////////////////////////////////////////
            for (std::size_t b = 0; b < triBlocks.size(); b++) {
                SPtr<Block3D> block = triBlocks[b];

                ////////////////////////////////////////////////////////////////////////////
                //// Block Dreieck-/test
                ////////////////////////////////////////////////////////////////////////////
                UbTupleDouble3 coords = grid.lock()->getBlockWorldCoordinates(block);
                UbTupleDouble3 deltas = grid.lock()->getBlockLengths(block);

                blockMinX[0] = (float)(val<1>(coords) - deltaMinX1[level]);
                blockMinX[1] = (float)(val<2>(coords) - deltaMinX2[level]);
                blockMinX[2] = (float)(val<3>(coords) - deltaMinX3[level]);

                blockMaxX[0] = (float)(val<1>(coords) + val<1>(deltas) + deltaMaxX1[level]);
                blockMaxX[1] = (float)(val<2>(coords) + val<2>(deltas) + deltaMaxX2[level]);
                blockMaxX[2] = (float)(val<3>(coords) + val<3>(deltas) + deltaMaxX3[level]);

                boxCenter[0] = (float)(0.5 * (blockMaxX[0] + blockMinX[0]));
                boxCenter[1] = (float)(0.5 * (blockMaxX[1] + blockMinX[1]));
                boxCenter[2] = (float)(0.5 * (blockMaxX[2] + blockMinX[2]));

                halfBoxSize[0] = (float)(0.5 * (blockMaxX[0] - blockMinX[0]));
                halfBoxSize[1] = (float)(0.5 * (blockMaxX[1] - blockMinX[1]));
                halfBoxSize[2] = (float)(0.5 * (blockMaxX[2] - blockMinX[2]));

                // wenn dreieck "vergroesserten cube" nicht schneidet/beruehrt -> keine BC moeglich -> continue
                if (!GbMeshTools3D::triBoxOverlap(boxCenter, halfBoxSize, triPoints)) {
                    counterTriBoxOverlap++;
                    continue;
                }

                //////////////////////////////////////////////////////////////////////////
                // Untersuchung der einzelnen nodes
                //////////////////////////////////////////////////////////////////////////
                bool blockGotBCs = false;

                SPtr<ILBMKernel> kernel  = block->getKernel();
                SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

                int indexMinX1 = 0;
                int indexMinX2 = 0;
                int indexMinX3 = 0;

                int indexMaxX1 = (int)bcMatrix->getNX1();
                int indexMaxX2 = (int)bcMatrix->getNX2();
                int indexMaxX3 = (int)bcMatrix->getNX3();

                std::set<std::vector<int>> &bcNodeIndices = this->bcNodeIndicesMap[block];
                //            std::set< UbTupleInt3 >& solidsFromOtherInteractors =
                //            tmpSolidNodesFromOtherInteractors[block];
                double q, distance;

                double &nodeDx1 = nodeDeltaToNeigh[level][DIR_P00];
                double &nodeDx2 = nodeDeltaToNeigh[level][DIR_0P0];
                double &nodeDx3 = nodeDeltaToNeigh[level][DIR_00P];

                // fuer OBB-Test
                double qEinflussDelta = 1.1 * sqrt(nodeDx1 * nodeDx1 + nodeDx2 * nodeDx2 + nodeDx3 * nodeDx3);

                for (int ix3 = indexMinX3; ix3 < indexMaxX3; ix3++) {
                    for (int ix2 = indexMinX2; ix2 < indexMaxX2; ix2++) {
                        for (int ix1 = indexMinX1; ix1 < indexMaxX1; ix1++) {
                            Vector3D pointplane1 = grid.lock()->getNodeCoordinates(block, ix1, ix2, ix3);
                            double internX1      = pointplane1[0];
                            double internX2      = pointplane1[1];
                            double internX3      = pointplane1[2];

                            //                     int blx1 = block->getX1();
                            //                     int blx2 = block->getX2();
                            //                     int blx3 = block->getX3();

                            if (bcMatrix->isSolid(ix1, ix2, ix3) || bcMatrix->isUndefined(ix1, ix2, ix3)) {
                                continue;
                            }

                            //////////////////////////////////////////////////////////////////////////
                            // Punkt in AABB von Dreieck?
                            //////////////////////////////////////////////////////////////////////////
                            // ehsan changed
                            bool pointIsOnBoundary = true;
                            if (!boundingCubeTriangle.isPointInGbObject3D(internX1, internX2, internX3,
                                                                          pointIsOnBoundary)) {
                                counterAABBTriFace++;
                                continue;
                            }
                            // std::cout<<"internX3  "<<internX3<<"  internX2"<<internX2<<" internX1 "<<internX1<<"\n";
                            //////////////////////////////////////////////////////////////////////////
                            // Halbebenentests
                            //////////////////////////////////////////////////////////////////////////
                            distance = halfSpace.getDistance(internX1, internX2, internX3);
                            // Punkt in Halbebene? (nein, wenn distance<0)
                            if (useHalfSpace &&
                                UbMath::less(distance, 0.0)) //== !halfSpace.ptInside(internX1,internX2,internX3) )
                            {
                                counterHalfspace++;
                                continue;
                            }

                            // BilligOBB-Test: wenn distance > qEinflussDelta -> kein q
                            if (UbMath::greater(fabs(distance), qEinflussDelta)) {
                                counterBilligOBB++;
                                continue;
                            }

                            /////////////////////////////////////////////////////////////////////////////
                            // Raytracingfür diskrete Boltzmannrichtungen
                            /////////////////////////////////////////////////////////////////////////////
                            gotQs = false;
                            bc    = SPtr<BoundaryConditions>();

                            // RAYTRACING - diskrete LB-dir zu Dreick
                            // e1 = v1 - v0
                            e1x1 = v2.x - v1.x;
                            e1x2 = v2.y - v1.y;
                            e1x3 = v2.z - v1.z;

                            // e2 = v2 - v0
                            e2x1 = v3.x - v1.x;
                            e2x2 = v3.y - v1.y;
                            e2x3 = v3.z - v1.z;

                            // s = o - v0
                            sx1 = internX1 - v1.x;
                            sx2 = internX2 - v1.y;
                            sx3 = internX3 - v1.z;

                            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                // p = d x e2
                                px1 = this->rayX2[fdir] * e2x3 - this->rayX3[fdir] * e2x2;
                                px2 = this->rayX3[fdir] * e2x1 - this->rayX1[fdir] * e2x3;
                                px3 = this->rayX1[fdir] * e2x2 - this->rayX2[fdir] * e2x1;

                                // a = e1 dot p
                                a = e1x1 * px1 + e1x2 * px2 + e1x3 * px3;
                                if (fabs(a) < 1.E-10)
                                    continue;
                                f = 1.0 / a;

                                // u = f * ( s dot p)
                                u = f * (sx1 * px1 + sx2 * px2 + sx3 * px3);
                                if (u < -1.E-10 || u > 1.0 + 1.E-10)
                                    continue;

                                // q = s x e1
                                qx1 = sx2 * e1x3 - sx3 * e1x2;
                                qx2 = sx3 * e1x1 - sx1 * e1x3;
                                qx3 = sx1 * e1x2 - sx2 * e1x1;

                                // v = f*(e2 dot q)
                                v = f * (this->rayX1[fdir] * qx1 + this->rayX2[fdir] * qx2 + this->rayX3[fdir] * qx3);
                                if (v < -1.E-10 || (u + v) > 1.0 + 1.E-10)
                                    continue;

                                // t = f * (e2 dot q)
                                q = f * (e2x1 * qx1 + e2x2 * qx2 + e2x3 * qx3);
                                q /= nodeDeltaToNeigh[level][fdir];
                                /////ehsan q/////////////////////////////////////////////////////////////////////
                                double det = triangle.nx * this->rayX1[fdir] + triangle.ny * this->rayX2[fdir] +
                                             triangle.nz * this->rayX3[fdir];

                                if (det > -1.E-10)
                                    continue;
                                double d  = triangle.nx * v1.x + triangle.ny * v1.y + triangle.nz * v1.z;
                                double x1 = -((-d * this->rayX1[fdir] - triangle.ny * this->rayX2[fdir] * internX1 -
                                               triangle.nz * this->rayX3[fdir] * internX1 +
                                               triangle.ny * this->rayX1[fdir] * internX2 +
                                               triangle.nz * this->rayX1[fdir] * internX3)) /
                                            det;
                                double y1 = -((-d * this->rayX2[fdir] + triangle.nx * this->rayX2[fdir] * internX1 -
                                               triangle.nx * this->rayX1[fdir] * internX2 -
                                               triangle.nz * this->rayX3[fdir] * internX2 +
                                               triangle.nz * this->rayX2[fdir] * internX3)) /
                                            det;
                                double z1 = -((-d * this->rayX3[fdir] + triangle.nx * this->rayX3[fdir] * internX1 +
                                               triangle.ny * this->rayX3[fdir] * internX2 -
                                               triangle.nx * this->rayX1[fdir] * internX3 -
                                               triangle.ny * this->rayX2[fdir] * internX3)) /
                                            det;
                                double q_ehsan =
                                    sqrt((x1 - internX1) * (x1 - internX1) + (y1 - internX2) * (y1 - internX2) +
                                         (z1 - internX3) * (z1 - internX3));
                                q_ehsan /= nodeDeltaToNeigh[level][fdir];
                                q = q_ehsan;
                                if (UbMath::greater(q, 1.0) || UbMath::lessEqual(q, 0.0))
                                    continue;

                                // gefundenes q auf gueltigkeit pruefen
                                if (UbMath::zero(q)) {
                                    // neu (18.05.2010)
                                    // es kann vorkommen, dass bei dünnwandigen geos punkte, die auf einem dreieck
                                    // liegen, qs bekommen, die durch die geo durchgehen. diese punkte werden später
                                    // jedoch nicht mehr auf solid getestet, da sie ja ne BC bekommen haben
                                    //--> da mind ein q==0.0 für eines der dreiecke -> dort solid setzen
                                    this->solidNodeIndicesMap[block].insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcMatrix->setSolid(ix1, ix2, ix3);
                                    continue;
                                }

                                if (UbMath::inClosedInterval(q, 1.0, 1.0))
                                    q = 1.0;
                                if (UbMath::greater(q, 0.0) && UbMath::lessEqual(q, 1.0)) {
                                    gotQs = blockGotBCs = true;

                                    bc = bcMatrix->getBC(ix1, ix2, ix3);

                                    // SG 26.08.2010 if(!bc && !bcMatrix->isSolid())
                                    if (!bc) {
                                        bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                                        ;
                                        bcMatrix->setBC(ix1, ix2, ix3, bc);
                                    } else if (UbMath::less(bc->getQ(fdir), q) &&
                                               UbMath::equal(-999.0, q)) // schon ein kuerzeres q voehanden?
                                    {
                                        // neu:: 18.05.2010
                                        // um falsche qs die evtl durch die "wand" gehen zu vermeiden
                                        // q nur dann neu setzen, wenn neues q kleiner als vorhandenes!
                                        // Zudem: insbesondere an ecken mit zwei BC geos ist nur das
                                        // naehere gueltig
                                        continue;
                                    }

                                    bc->setBoundaryVelocityX1(vx1);
                                    bc->setBoundaryVelocityX2(vx2);
                                    bc->setBoundaryVelocityX3(vx3);

                                    for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                        this->BCs[index]->adaptBCForDirection(*this, bc, internX1, internX2,
                                                                                     internX3, q, fdir);

                                    // fuer beschleunigtes wiedereinlesen
                                    if (this->reinitWithStoredQsFlag) {
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)].resize(
                                            D3Q27System::FENDDIR + 1 + 3, -1.0f);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)][fdir] = float(q);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 0] = float(internX1);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 1] = float(internX2);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 2] = float(internX3);
                                    }
                                }
                            }

                            if (gotQs) {
                                std::vector<int> p(3);
                                p[0] = ix1;
                                p[1] = ix2;
                                p[2] = ix3;
                                bcNodeIndices.insert(p);

                                for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                    this->BCs[index]->adaptBC(*this, bc, internX1, internX2, internX3);
                            }
                        }
                    }
                }
            }
            // dynamische Punkte des GbCuboids muessen leider per "Hand" geloescht werden :-(
            boundingCubeTriangle.finalize();
        }
    }
    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - setQs end ");
}
//////////////////////////////////////////////////////////////////////////
// Vorgehesnweise
// A – Bestimmung der q's
//  1. fuer jeden Bounding cube eines Dreiecks des netzes werden die Bloecke des Blockgitter ermittelt
//  2. mittels eines Dreieck/Block Verschneidungstest werden weitere nicht relevante Bloecke aussortiert
//     (fuer lange „schief“ im Raum stehende Dreicke, bei denen das Bounding Cube suboptimal ist)
//  3. jeder Knoten dieser blöcke wird gegen das bound cube des dreiecks getestet
//  4. Knoten die innerhalb des Cubes aber „innerhalb“ des Netzes liegen werden mittels Halbebenentest aussoriert
//  5. fuer die restliche Knoten erfolgt die q bestimmung mittels effizienter raytracing algorithmen
//     fuer die diskreten Boltzmannrichtungen
// B – Setzen der nicht aktiven Bloecke und Solid Nodes
//  alle Bloecke des Bounding Cube des Netzes, die mind eine BC erhielten, wurden in A markiert
//  1. fuer nicht markierte Bloecke genuegt EIN pointInObject(Dreicksnetz)-Test um den gesamten Block bei Erfolg als
//  „not active“ zu markieren
//  2. fuer markiertre Bloecke wird ein rekursiver Fuellalgorithmus durchgefuehrt
void D3Q27TriFaceMeshInteractor::initInteractor2(const real &timeStep)
{
    using namespace vf::lbm::dir;

    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - initInteractor start ");
    if (!this->grid.lock())
        throw UbException(UB_EXARGS, "ups, no grid.lock()!!");

    if (this->reinitWithStoredQsFlag && !bcNodeIndicesAndQsMap.empty()) {
        this->reinitWithStoredQs(timeStep);
        return;
    }

    GbTriFaceMesh3D *mesh = dynamic_cast<GbTriFaceMesh3D *>(this->geoObject3D.get());

    UBLOGML(logDEBUG1, "\nLBMTriFaceMeshInteractor - initInteractor for \"" << mesh->getName() << " \" t=" << timeStep);
    // cout<<" - init basics ...";

    this->removeBcBlocks(); // hier wird auch die nodeIndicesMap geloescht!
    this->removeSolidBlocks();

    //////////////////////////////////////////////////////////////////////////
    // init bcs
    //////////////////////////////////////////////////////////////////////////
    int nofAdapter = (int)this->BCs.size();
    if (nofAdapter == 0)
        std::cout
            << "WARNING - D3Q27TriFaceMeshInteractor::initInteractor Warning - no nodeAdapter available for " /*<<this->getName()*/
            << std::endl;
    bool needTimeDependence = false;
    for (int pos = 0; pos < nofAdapter; ++pos) {
        this->BCs[pos]->init(this, timeStep);
        if (this->BCs[pos]->isTimeDependent())
            needTimeDependence = true;
    }
    if (needTimeDependence)
        this->setTimeDependent();
    else
        this->unsetTimeDependent();

    //////////////////////////////////////////////////////////////////////////
    // grid.lock() info
    //////////////////////////////////////////////////////////////////////////
    int coarsestInitLevel = grid.lock()->getCoarsestInitializedLevel();
    int finestInitLevel   = grid.lock()->getFinestInitializedLevel();

    UbTupleInt3 blocknx = grid.lock()->getBlockNX();
    int blocknx1        = val<1>(blocknx); // gilt fuer alle Level
    int blocknx2        = val<2>(blocknx); // gilt fuer alle Level
    int blocknx3        = val<3>(blocknx); // gilt fuer alle Level

    // grobe Blocklaengen
    SPtr<CoordinateTransformation3D> trafo = grid.lock()->getCoordinateTransformator();
    double cblockDeltaX1, cblockDeltaX2, cblockDeltaX3, delta;
    cblockDeltaX1 = cblockDeltaX2 = cblockDeltaX3 = delta = 1.0 / (double)(1 << coarsestInitLevel);
    if (trafo) {
        cblockDeltaX1 = trafo->getX1CoordinateScaling() * delta;
        cblockDeltaX2 = trafo->getX2CoordinateScaling() * delta;
        cblockDeltaX3 = trafo->getX3CoordinateScaling() * delta;
    }
    // levelspezifische blocklaengen und knotenabstaende
    std::vector<std::vector<double>> nodeDeltaToNeigh(finestInitLevel + 1);
    // vector<double> blockDeltaX1(finestInitLevel+1), blockDeltaX2(finestInitLevel+1), blockDeltaX3(finestInitLevel+1);
    std::vector<float> deltaMinX1(finestInitLevel + 1), deltaMinX2(finestInitLevel + 1),
        deltaMinX3(finestInitLevel + 1);
    std::vector<float> deltaMaxX1(finestInitLevel + 1), deltaMaxX2(finestInitLevel + 1),
        deltaMaxX3(finestInitLevel + 1);

    // Im Boltzmankontext muss dx1==dx2==dx3 sein!!
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX2 / (double)blocknx2));
    assert(UbMath::equal(cblockDeltaX1 / (double)blocknx1, cblockDeltaX3 / (double)blocknx3));

    for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
        double nodeDelta = cblockDeltaX1 / (double)(blocknx1 * (1 << (level - coarsestInitLevel)));

        std::vector<double> distNeigh(D3Q27System::FENDDIR + 1, 0.0);
        D3Q27System::calcDistanceToNeighbors(distNeigh, nodeDelta);

        nodeDeltaToNeigh[level].resize(D3Q27System::ENDDIR + 1, 0.0);
        for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
            nodeDeltaToNeigh[level][fdir] = distNeigh[fdir];
        }

        // im gegensatz zum allg. Cell3DInteractor kann man hier auf max(0.02*blockDeltaX1[level],fabs(...)) verzichten
        // da dies nur für blockDeltaCalculator->getMinX1Delta(level)==0.0 benötigt wird. ist im D3Q19... aber nie so
        // Geller: kann man nicht diesen befuckten DeltaCalculator weglassen und hier einfach die Formel zum Delta
        // rechnen reinpacken SirAnn: klar, mann kann auch weißwuerste am Alex verkaufen... aber zum einen ist das Ding
        // dazu da
        // und zum anderen sollt eman mal überlegen: "Formel zum Delta rechnen"->man muss rechnen
        // blockDeltaCalculator->getMinX1Delta(level) -> ein geinlinter wert wird geholt

        deltaMinX1[level] = (float)(1.2 * nodeDelta); // kein minus da unten -deltaMin
        deltaMinX2[level] = (float)(1.2 * nodeDelta);
        deltaMinX3[level] = (float)(1.2 * nodeDelta);
        deltaMaxX1[level] = (float)(1.2 * nodeDelta);
        deltaMaxX2[level] = (float)(1.2 * nodeDelta);
        deltaMaxX3[level] = (float)(1.2 * nodeDelta);
    }

    //////////////////////////////////////////////////////////////////////////
    // bounding cubes des TriFaceMesh ermitteln (pro level)
    //////////////////////////////////////////////////////////////////////////
    // min/max Werte des Dreiecksnetzes holen
    double geoMinX1(0.0), geoMinX2(0.0), geoMinX3(0.0), geoMaxX1(0.0), geoMaxX2(0.0), geoMaxX3(0.0);
    if (this->isSolid() || this->isMoveable()) {
        geoMinX1 = this->geoObject3D->getX1Minimum();
        geoMaxX1 = this->geoObject3D->getX1Maximum();
        geoMinX2 = this->geoObject3D->getX2Minimum();
        geoMaxX2 = this->geoObject3D->getX2Maximum();
        geoMinX3 = this->geoObject3D->getX3Minimum();
        geoMaxX3 = this->geoObject3D->getX3Maximum();
    } else
        throw UbException(UB_EXARGS, "only TYPE==SOLID is implemented");

    std::map<SPtr<Block3D>, SolidCheckMethod> blocksForSolidCheck;

    for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
        if (this->isSolid() || this->isMoveable()) {
            // bloecke fuer "bounding cube gesamt"
            std::vector<SPtr<Block3D>> tmpblocks;
            grid.lock()->getBlocksByCuboid(level, geoMinX1 - deltaMinX1[level], geoMinX2 - deltaMinX2[level],
                                           geoMinX3 - deltaMinX3[level], geoMaxX1 + deltaMaxX1[level],
                                           geoMaxX2 + deltaMaxX2[level], geoMaxX3 + deltaMaxX3[level], tmpblocks);

            for (size_t i = 0; i < tmpblocks.size(); i++)
                blocksForSolidCheck[tmpblocks[i]] = PointInObject;
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // FE-specific
    //////////////////////////////////////////////////////////////////////////
    // bool calcVelocities = false;
    // FeTriFaceMesh3D* feMesh = dynamic_cast<FeTriFaceMesh3D*>(mesh);
    // std::vector<FeTriFaceMesh3D::VertexAttributes>* attributes = NULL;
    // if(feMesh)
    //{
    //   calcVelocities = true;
    //   attributes     = feMesh->getAttributes();
    //}

    //////////////////////////////////////////////////////////////////////////
    // DREIECKE: q-Bestimmung
    //////////////////////////////////////////////////////////////////////////

    // notwendige variablen initialisieren (u.a. blockDeltas des groben levels)
    float triPoints[3][3];
    real vx1 = 0.0, vx2 = 0.0, vx3 = 0.0;
    unsigned counterTriBoxOverlap = 0, counterAABBTriFace = 0, counterHalfspace = 0, counterBilligOBB = 0;
    std::vector<GbTriFaceMesh3D::TriFace> &triangles = *mesh->getTriangles();
    std::vector<GbTriFaceMesh3D::Vertex> &nodes      = *mesh->getNodes();
    std::map<SPtr<Block3D>, std::set<std::vector<int>>> tmpSolidNodesFromOtherInteractors;

    int onePercent = UbMath::integerRounding(triangles.size() * 0.01);
    if (onePercent == 0)
        onePercent = 1;
    UbTimer setQTimer;
    setQTimer.start();
    UBLOG(logDEBUG3, " - setQs for " << (int)triangles.size() << " triangles");

    //   bool solidFromOtherInteractor = false;
    float blockMinX[3], blockMaxX[3], boxCenter[3], halfBoxSize[3];

    for (size_t t = 0; t < triangles.size(); t++) {
        // if (t==10577)
        //{
        // int ehsan=0;
        //}
        //////////////////////////////////////////////////////////////////////////
        // Halfspace zum Dreieck generieren und min/max des Dreiecks ermitteln
        //////////////////////////////////////////////////////////////////////////
        GbTriFaceMesh3D::TriFace &triangle = triangles[t];

        GbTriFaceMesh3D::Vertex &v1 = nodes[triangle.v1];
        GbTriFaceMesh3D::Vertex &v2 = nodes[triangle.v2];
        GbTriFaceMesh3D::Vertex &v3 = nodes[triangle.v3];

        GbHalfSpace3D halfSpace(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z);

        // if(calcVelocities)
        //{
        //   FeTriFaceMesh3D::VertexAttributes& vAttribut1 = (*attributes)[triangle.v1];
        //   FeTriFaceMesh3D::VertexAttributes& vAttribut2 = (*attributes)[triangle.v2];
        //   FeTriFaceMesh3D::VertexAttributes& vAttribut3 = (*attributes)[triangle.v3];
        //   vx1 =
        //   (float)(UbMath::c1o3*(vAttribut1.getVelocityX()+vAttribut2.getVelocityX()+vAttribut3.getVelocityX())); vx2
        //   = (float)(UbMath::c1o3*(vAttribut1.getVelocityY()+vAttribut2.getVelocityY()+vAttribut3.getVelocityY()));
        //   vx3 =
        //   (float)(UbMath::c1o3*(vAttribut1.getVelocityZ()+vAttribut2.getVelocityZ()+vAttribut3.getVelocityZ()));
        //}

        //////////////////////////////////////////////////////////////////////////
        // fuer GbMeshTools3D::triBoxOverlap
        //////////////////////////////////////////////////////////////////////////
        triPoints[0][0] = v1.x;
        triPoints[0][1] = v1.y;
        triPoints[0][2] = v1.z;
        triPoints[1][0] = v2.x;
        triPoints[1][1] = v2.y;
        triPoints[1][2] = v2.z;
        triPoints[2][0] = v3.x;
        triPoints[2][1] = v3.y;
        triPoints[2][2] = v3.z;

        double minX1 = triangle.getMinX(nodes);
        double maxX1 = triangle.getMaxX(nodes);
        double minX2 = triangle.getMinY(nodes);
        double maxX2 = triangle.getMaxY(nodes);
        double minX3 = triangle.getMinZ(nodes);
        double maxX3 = triangle.getMaxZ(nodes);

        //////////////////////////////////////////////////////////////////////////
        // Schleife ueber alle Level
        //////////////////////////////////////////////////////////////////////////
        double e1x1, e1x2, e1x3, e2x1, e2x2, e2x3, px1, px2, px3, a, f, sx1, sx2, sx3, u, qx1, qx2, qx3, v;
        bool gotQs = false;
        SPtr<BoundaryConditions> bc;

        for (int level = coarsestInitLevel; level <= finestInitLevel; level++) {
            //////////////////////////////////////////////////////////////////////////
            // levelspezifisches BoundCube des Dreicks ermitteln und zugehörige Bloecke beziehen
            //////////////////////////////////////////////////////////////////////////
            double boundCubeTriangleMinX1 = minX1 - deltaMinX1[level];
            double boundCubeTriangleMaxX1 = maxX1 + deltaMaxX1[level];
            double boundCubeTriangleMinX2 = minX2 - deltaMinX2[level];
            double boundCubeTriangleMaxX2 = maxX2 + deltaMaxX2[level];
            double boundCubeTriangleMinX3 = minX3 - deltaMinX3[level];
            double boundCubeTriangleMaxX3 = maxX3 + deltaMaxX3[level];

            GbCuboid3D boundingCubeTriangle(boundCubeTriangleMinX1, boundCubeTriangleMinX2, boundCubeTriangleMinX3,
                                            boundCubeTriangleMaxX1, boundCubeTriangleMaxX2, boundCubeTriangleMaxX3);

            std::vector<SPtr<Block3D>> triBlocks;
            grid.lock()->getBlocksByCuboid(level, boundCubeTriangleMinX1, boundCubeTriangleMinX2,
                                           boundCubeTriangleMinX3, boundCubeTriangleMaxX1, boundCubeTriangleMaxX2,
                                           boundCubeTriangleMaxX3, triBlocks);

            //////////////////////////////////////////////////////////////////////////
            // Schleife ueber bloecke des level, die das dreieck beinhalten
            //////////////////////////////////////////////////////////////////////////
            for (std::size_t b = 0; b < triBlocks.size(); b++) {
                SPtr<Block3D> block = triBlocks[b];

                ////////////////////////////////////////////////////////////////////////////
                //// Block Dreieck-/test
                ////////////////////////////////////////////////////////////////////////////
                UbTupleDouble3 coords = grid.lock()->getBlockWorldCoordinates(block);
                UbTupleDouble3 deltas = grid.lock()->getBlockLengths(block);

                blockMinX[0] = (float)(val<1>(coords) - deltaMinX1[level]);
                blockMinX[1] = (float)(val<2>(coords) - deltaMinX2[level]);
                blockMinX[2] = (float)(val<3>(coords) - deltaMinX3[level]);

                blockMaxX[0] = (float)(val<1>(coords) + val<1>(deltas) + deltaMaxX1[level]);
                blockMaxX[1] = (float)(val<2>(coords) + val<2>(deltas) + deltaMaxX2[level]);
                blockMaxX[2] = (float)(val<3>(coords) + val<3>(deltas) + deltaMaxX3[level]);

                boxCenter[0] = (float)(0.5 * (blockMaxX[0] + blockMinX[0]));
                boxCenter[1] = (float)(0.5 * (blockMaxX[1] + blockMinX[1]));
                boxCenter[2] = (float)(0.5 * (blockMaxX[2] + blockMinX[2]));

                halfBoxSize[0] = (float)(0.5 * (blockMaxX[0] - blockMinX[0]));
                halfBoxSize[1] = (float)(0.5 * (blockMaxX[1] - blockMinX[1]));
                halfBoxSize[2] = (float)(0.5 * (blockMaxX[2] - blockMinX[2]));

                // wenn dreieck "vergroesserten cube" nicht schneidet/beruehrt -> keine BC moeglich -> continue
                if (!GbMeshTools3D::triBoxOverlap(boxCenter, halfBoxSize, triPoints)) {
                    counterTriBoxOverlap++;
                    continue;
                }

                //////////////////////////////////////////////////////////////////////////
                // Untersuchung der einzelnen nodes
                //////////////////////////////////////////////////////////////////////////
                bool blockGotBCs = false;

                SPtr<ILBMKernel> kernel  = block->getKernel();
                SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

                int indexMinX1 = 0;
                int indexMinX2 = 0;
                int indexMinX3 = 0;

                int indexMaxX1 = (int)bcMatrix->getNX1();
                int indexMaxX2 = (int)bcMatrix->getNX2();
                int indexMaxX3 = (int)bcMatrix->getNX3();

                std::set<std::vector<int>> &bcNodeIndices              = this->bcNodeIndicesMap[block];
                std::set<std::vector<int>> &solidsFromOtherInteractors = tmpSolidNodesFromOtherInteractors[block];
                double q, internX1, internX2, internX3, distance;

                double &nodeDx1 = nodeDeltaToNeigh[level][DIR_P00];
                double &nodeDx2 = nodeDeltaToNeigh[level][DIR_0P0];
                double &nodeDx3 = nodeDeltaToNeigh[level][DIR_00P];

                // fuer OBB-Test
                double qEinflussDelta = 1.1 * sqrt(nodeDx1 * nodeDx1 + nodeDx2 * nodeDx2 + nodeDx3 * nodeDx3);

                for (int ix3 = indexMinX3; ix3 < indexMaxX3; ix3++) {
                    internX3 = val<3>(coords) + nodeDx3 * ix3 - 0.5 * nodeDx3;
                    for (int ix2 = indexMinX2; ix2 < indexMaxX2; ix2++) {
                        internX2 = val<2>(coords) + nodeDx2 * ix2 - 0.5 * nodeDx2;
                        for (int ix1 = indexMinX1; ix1 < indexMaxX1; ix1++) {

                            //					  int blx1 =block->getX1();
                            //					  int blx2 = block->getX2();
                            //					  int blx3 = block->getX3();

                            //					  if (blx1==0&&blx2==1&&blx3==0)
                            //					  {
                            //						  //if (ix2==39&&ix3==4)
                            //							   if (ix2==39&&ix3==4)
                            //						  {
                            //							 int seb=0;
                            //						  }
                            //					  }
                            // Problem: wenn voher der punkt durch eine andere geo not active gesetzt wird und
                            // dieser nun uebersprungen wird, dann hat man spaeter beim fuellalgorithmus luecken
                            // in der front und der block wird u.U. faelschlicher weise komplett solid markiert
                            // Lsg: positionen merken und erst Nach dem fuellarlgo wieder auf not active setzen :-)
                            //                     solidFromOtherInteractor = false;
                            if (bcMatrix->isSolid(ix1, ix2, ix3)) {
                                if (this->reinitWithStoredQsFlag) {
                                    // solidFromOtherInteractor = true;   //hier muss man weitermachen
                                    // SG //oje
                                    std::vector<int> p(3);
                                    p[0] = ix1;
                                    p[1] = ix2;
                                    p[2] = ix3;
                                    solidsFromOtherInteractors.insert(p);
                                } else {
                                    // SG //oje
                                    std::vector<int> p(3);
                                    p[0] = ix1;
                                    p[1] = ix2;
                                    p[2] = ix3;
                                    solidsFromOtherInteractors.insert(p);
                                    // SG continue;
                                    // solidFromOtherInteractor = true;
                                }
                            }

                            internX1 = val<1>(coords) + nodeDx1 * ix1 - 0.5 * nodeDx1;

                            //////////////////////////////////////////////////////////////////////////
                            // Punkt in AABB von Dreieck?
                            //////////////////////////////////////////////////////////////////////////
                            // ehsan changed◘
                            bool pointIsOnBoundary = false;
                            if (!boundingCubeTriangle.isPointInGbObject3D(internX1, internX2, internX3,
                                                                          pointIsOnBoundary))
                            // if( !boundingCubeTriangle.isPointInGbObject3D(internX1, internX2, internX3) )
                            {
                                counterAABBTriFace++;
                                continue;
                            }

                            //////////////////////////////////////////////////////////////////////////
                            // Halbebenentests
                            //////////////////////////////////////////////////////////////////////////
                            distance = halfSpace.getDistance(internX1, internX2, internX3);

                            // Punkt in Halbebene? (nein, wenn distance<0)
                            if (useHalfSpace &&
                                UbMath::less(distance, 0.0)) //== !halfSpace.ptInside(internX1,internX2,internX3) )
                            {
                                counterHalfspace++;
                                continue;
                            }

                            // BilligOBB-Test: wenn distance > qEinflussDelta -> kein q
                            if (UbMath::greater(fabs(distance), qEinflussDelta)) {
                                counterBilligOBB++;
                                continue;
                            }

                            /////////////////////////////////////////////////////////////////////////////
                            // Raytracingfür diskrete Boltzmannrichtungen
                            /////////////////////////////////////////////////////////////////////////////
                            gotQs = false;
                            bc    = SPtr<BoundaryConditions>();

                            // RAYTRACING - diskrete LB-dir zu Dreick
                            // e1 = v1 - v0
                            e1x1 = v2.x - v1.x;
                            e1x2 = v2.y - v1.y;
                            e1x3 = v2.z - v1.z;

                            // e2 = v2 - v0
                            e2x1 = v3.x - v1.x;
                            e2x2 = v3.y - v1.y;
                            e2x3 = v3.z - v1.z;

                            // s = o - v0
                            sx1 = internX1 - v1.x;
                            sx2 = internX2 - v1.y;
                            sx3 = internX3 - v1.z;

                            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                                // p = d x e2
                                px1 = this->rayX2[fdir] * e2x3 - this->rayX3[fdir] * e2x2;
                                px2 = this->rayX3[fdir] * e2x1 - this->rayX1[fdir] * e2x3;
                                px3 = this->rayX1[fdir] * e2x2 - this->rayX2[fdir] * e2x1;

                                // a = e1 dot p
                                a = e1x1 * px1 + e1x2 * px2 + e1x3 * px3;
                                if (fabs(a) < 1.E-10)
                                    continue;
                                f = 1.0 / a;

                                // u = f * ( s dot p)
                                u = f * (sx1 * px1 + sx2 * px2 + sx3 * px3);
                                if (u < -1.E-10 || u > 1.0 + 1.E-10)
                                    continue;

                                // q = s x e1
                                qx1 = sx2 * e1x3 - sx3 * e1x2;
                                qx2 = sx3 * e1x1 - sx1 * e1x3;
                                qx3 = sx1 * e1x2 - sx2 * e1x1;

                                // v = f*(e2 dot q)
                                v = f * (this->rayX1[fdir] * qx1 + this->rayX2[fdir] * qx2 + this->rayX3[fdir] * qx3);
                                if (v < -1.E-10 || (u + v) > 1.0 + 1.E-10)
                                    continue;

                                // t = f * (e2 dot q)
                                q = f * (e2x1 * qx1 + e2x2 * qx2 + e2x3 * qx3);
                                q /= nodeDeltaToNeigh[level][fdir];

                                // gefundenes q auf gueltigkeit pruefen
                                if (UbMath::zero(q)) {
                                    // neu (18.05.2010)
                                    // es kann vorkommen, dass bei dünnwandigen geos punkte, die auf einem dreieck
                                    // liegen, qs bekommen, die durch die geo durchgehen. diese punkte werden später
                                    // jedoch nicht mehr auf solid getestet, da sie ja ne BC bekommen haben
                                    //--> da mind ein q==0.0 für eines der dreiecke -> dort solid setzen
                                    this->solidNodeIndicesMap[block].insert(UbTupleInt3(ix1, ix2, ix3));
                                    bcMatrix->setSolid(ix1, ix2, ix3);
                                    continue;
                                }

                                if (UbMath::inClosedInterval(q, 1.0, 1.0))
                                    q = 1.0;
                                if (UbMath::greater(q, 0.0) && UbMath::lessEqual(q, 1.0)) {
                                    // if( !solidFromOtherInteractor ) //--> Knoten schon solid-->BC setzen
                                    // ueberfluessig SG changed to if( solidFromOtherInteractor ) //--> Knoten schon
                                    // solid-->BC setzen ueberfluessig
                                    {
                                        // SG 26.08.2010 muss bereits hierhin, da das continue sonst den Knoten nicht
                                        // als transNode fürs markiert
                                        gotQs = blockGotBCs = true;

                                        bc = bcMatrix->getBC(ix1, ix2, ix3);

                                        // SG 26.08.2010 if(!bc && !bcMatrix->isSolid())
                                        if (!bc) {
                                            bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                                            ;
                                            bcMatrix->setBC(ix1, ix2, ix3, bc);
                                        } else if (UbMath::less(bc->getQ(fdir), q)) // schon ein kuerzeres q voehanden?
                                        {
                                            // neu:: 18.05.2010
                                            // um falsche qs die evtl durch die "wand" gehen zu vermeiden
                                            // q nur dann neu setzen, wenn neues q kleiner als vorhandenes!
                                            // Zudem: insbesondere an ecken mit zwei BC geos ist nur das
                                            // naehere gueltig
                                            continue;
                                        }

                                        bc->setBoundaryVelocityX1(vx1);
                                        bc->setBoundaryVelocityX2(vx2);
                                        bc->setBoundaryVelocityX3(vx3);

                                        for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                            this->BCs[index]->adaptBCForDirection(*this, bc, internX1, internX2,
                                                                                         internX3, q, fdir);

                                        // SG 26.08.2010 gotQs=blockGotBCs=true;
                                    }
                                    // fuer beschleunigtes wiedereinlesen
                                    if (this->reinitWithStoredQsFlag) {
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)].resize(
                                            D3Q27System::FENDDIR + 1 + 3, -1.0f);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)][fdir] = float(q);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 0] = float(internX1);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 1] = float(internX2);
                                        bcNodeIndicesAndQsMap[block][UbTupleInt3(ix1, ix2, ix3)]
                                                             [D3Q27System::FENDDIR + 1 + 2] = float(internX3);
                                    }
                                }
                            }

                            if (gotQs) {
                                std::vector<int> p(3);
                                p[0] = ix1;
                                p[1] = ix2;
                                p[2] = ix3;
                                bcNodeIndices.insert(p);

                                for (int index = (int)this->BCs.size() - 1; index >= 0; --index)
                                    this->BCs[index]->adaptBC(*this, bc, internX1, internX2, internX3);
                            }
                        }
                    }
                }
                // Block wird für scanline-check "anmelden", dieser wird dann spaeter zu "transBlocks" hinzugefuegt
                if (blockGotBCs) {
                    blocksForSolidCheck[block] = ScanLine;
                }
                //            bvd->getTimer().stop();
            }
            // dynamische Punkte des GbCuboids muessen leider per "Hand" geloescht werden :-(
            boundingCubeTriangle.finalize();
        }
    }
    setQTimer.stop();

    UBLOG(logDEBUG1, " - setQs for " /*<< this->getName()*/ << " - " << (int)triangles.size()
                                                            << " triangles: 100% done in " << setQTimer.getTotalTime()
                                                            << "sec");
    UBLOG(logDEBUG1, "       * rejected blocks with tribox overlap test : " << counterTriBoxOverlap);
    UBLOG(logDEBUG1, "       * rejected nodes  with AABB           test : " << counterAABBTriFace);
    UBLOG(logDEBUG1, "       * rejected nodes  with halfspace      test : " << counterHalfspace);
    UBLOG(logDEBUG1, "       * rejected nodes  with OBB            test : " << counterBilligOBB);

    typedef std::map<SPtr<Block3D>, SolidCheckMethod>::iterator BlockSolidCheckMethodIterator;

    //////////////////////////////////////////////////////////////////////////
    // SOLID checks
    //////////////////////////////////////////////////////////////////////////
    if (regardPIOTest) {
        int pointInObjectCounter = 0;
        int scanlineCounter      = 0;
        //      int counter              = 0;

        // sollte die matrix groesse zu gross sein und der rekursive floodFill mehr speicher
        // benoetigen als der Stack hergibt -> keinen floodFill verwenden!
        void (D3Q27TriFaceMeshInteractor::*gridFill)(CbArray3D<FLAGS> &, const short &, const short &, const short &,
                                                     const FLAGS &) = NULL;
        /*if(blocknx1*blocknx2*blocknx3 < 200000 ) gridFill = &D3Q27TriFaceMeshInteractor::recursiveGridFill;*/
        /*else */ gridFill = &D3Q27TriFaceMeshInteractor::iterativeGridFill;

        UBLOG(logDEBUG1, " - setSolids for " << blocksForSolidCheck.size() << " blocks");

        UbTimer scanLineTimer;
        UbTimer solidTimer;
        solidTimer.start();

        for (BlockSolidCheckMethodIterator pos = blocksForSolidCheck.begin(); pos != blocksForSolidCheck.end(); ++pos) {
            SPtr<Block3D> const &block = pos->first;
            int level                  = block->getLevel();

            UbTupleDouble3 coords = grid.lock()->getBlockWorldCoordinates(block);

            // Bloecke, die keinerlei Verschneidung mit Dreicken bzw. deren BoundCubes hatten
            // hier: durch inside/outside Tests EINES knotens gilt fuer ALLE Knoten des blockes
            if (pos->second == PointInObject) {
                pointInObjectCounter++;
                if (mesh->isPointInGbObject3D(val<1>(coords), val<2>(coords), val<3>(coords))) {
                    // block->setActive(false);
                    // this->solidBlocks.push_back(block);
                }
            }
            // Bloecke, die Verschneidung mit Dreicken bzw. deren BoundCubes hatten
            // scanline algortihmus. dieser berücksichtigt durch weitere tests, dass innerhalb evtl schon andere
            // geos bcs gesetzt haben. es werden ausschließelich solids gesetzt (also keine FLUIDS oder neuen BCs)
            else if (pos->second == ScanLine) {
                scanlineCounter++;
                scanLineTimer.start();

                SPtr<ILBMKernel> kernel = block->getKernel();
                if (!kernel)
                    throw UbException(UB_EXARGS, "na sowas kein kernel bzw. kernel=NULL (2)");
                SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

                //            bvd->getTimer().start();
                //            int indexMinX1 = 0;
                //            int indexMinX2 = 0;
                //            int indexMinX3 = 0;
                int indexMaxX1 = (int)bcMatrix->getNX1();
                int indexMaxX2 = (int)bcMatrix->getNX2();
                int indexMaxX3 = (int)bcMatrix->getNX3();

                // quick and dirty
                blocknx1 = indexMaxX1;
                blocknx2 = indexMaxX2;
                blocknx3 = indexMaxX3;

                std::set<UbTupleInt3> &solidNodeIndices = this->solidNodeIndicesMap[block];

                float nodeDeltaX1 = (float)nodeDeltaToNeigh[level][DIR_P00];
                float nodeDeltaX2 = (float)nodeDeltaToNeigh[level][DIR_0P0];
                float nodeDeltaX3 = (float)nodeDeltaToNeigh[level][DIR_00P];

                // flagfield matrix initialisieren
                CbArray3D<FLAGS> flagField(blocknx1, blocknx2, blocknx3, UNDEF_FLAG);

                // hier gesetzte bcs markieren
                std::set<std::vector<int>> &transNodeIndices = this->bcNodeIndicesMap[block];
                std::set<std::vector<int>>::iterator setPos;
                for (setPos = transNodeIndices.begin(); setPos != transNodeIndices.end(); ++setPos)
                    flagField((*setPos)[0], (*setPos)[1], (*setPos)[2]) = BC_FLAG;

                // solids die bereits durch andere interaktoren gesetzt wurden (wurden oben gespeichert)
                // ist EMPTY bei reinitWithStoredQsFlag == true
                // SG 28.08.2010            std::set< UbTupleInt3 >& tmpSolidNodeIndices =
                // tmpSolidNodesFromOtherInteractors[block]; SG 28.08.2010  if(reinitWithStoredQsFlag &&
                // !tmpSolidNodeIndices.empty() ) throw UbException(UB_EXARGS, "tmpSolidNodeIndices darf bei
                // reinitWithStoredQsFlag==true keine Knoten enthalten"); SG 28.08.2010
                // for(setPos=tmpSolidNodeIndices.begin(); setPos!=tmpSolidNodeIndices.end();  ++setPos) SG 28.08.2010
                // flagField( val<1>(*setPos), val<2>(*setPos), val<3>(*setPos) ) = OLDSOLID_FLAG;

                // flagfield matrix belegen
                for (int bx3 = 0; bx3 < blocknx3; ++bx3) {
                    for (int bx2 = 0; bx2 < blocknx2; ++bx2) {
                        for (int bx1 = 0; bx1 < blocknx1; ++bx1) {

                            //					  if (bx2==9&&bx3==29)
                            //					  {
                            //						  int ride=0;
                            //					  }
                            if (flagField(bx1, bx2, bx3) == UNDEF_FLAG) {
                                if (mesh->isPointInGbObject3D(val<1>(coords) + bx1 * nodeDeltaX1 - 0.5 * nodeDeltaX1,
                                                              val<2>(coords) + bx2 * nodeDeltaX2 - 0.5 * nodeDeltaX2,
                                                              val<3>(coords) + bx3 * nodeDeltaX3 - 0.5 * nodeDeltaX3)) {
                                    (this->*gridFill)(flagField, bx1, bx2, bx3, SOLID_FLAG);
                                } else {
                                    (this->*gridFill)(flagField, bx1, bx2, bx3, FLUID_FLAG);
                                }
                            }

                            if (flagField(bx1, bx2, bx3) == SOLID_FLAG) {
                                // hier ist noch das Problem, das "alle" solid in die solidNodeIndices kommen
                                // evtl. Abhilfe durch einführen eines anderen Flags bei tmpSolidNodeIndices ..
                                solidNodeIndices.insert(UbTupleInt3(bx1, bx2, bx3));
                                bcMatrix->setSolid(bx1, bx2, bx3);
                            }
                            // SG 28.08.2010  else if( flagField(bx1,bx2,bx3)==OLDSOLID_FLAG )
                            // SG 28.08.2010  {
                            // SG 28.08.2010     bcMatrix->setSolid(bx1,bx2,bx3);
                            // SG 28.08.2010  }
                        }
                    }
                }

                // SG 28.08.2010 halt danach setzen, damit die BCs die fälschlicherweise gesetzt wurden korrigiert
                // werden
                std::set<std::vector<int>> &tmpSolidNodeIndices = tmpSolidNodesFromOtherInteractors[block];
                for (setPos = tmpSolidNodeIndices.begin(); setPos != tmpSolidNodeIndices.end(); ++setPos)
                    bcMatrix->setSolid((*setPos)[0], (*setPos)[1], (*setPos)[2]);

                // block hat  in initInteractor mind eine BC erhalten -> transBlock
                this->bcBlocks.push_back(block);
                scanLineTimer.stop();

                //            bvd->getTimer().stop();
            } else
                throw UbException(UB_EXARGS, "unknown option for in object test");
        }

        solidTimer.stop();

        UBLOG(logDEBUG1, " - setSolids for " << blocksForSolidCheck.size() << " blocks: 100% done in "
                                             << solidTimer.getTotalTime() << "s");
        UBLOG(logDEBUG1, "       * pointInObject for " << pointInObjectCounter << " blocks in "
                                                       << solidTimer.getTotalTime() - scanLineTimer.getTotalTime()
                                                       << "s");
        UBLOG(logDEBUG1, "       * flood fill    for " << scanlineCounter << " blocks in "
                                                       << scanLineTimer.getTotalTime() << " secs");
        UBLOG(logDEBUG1, "LBMTriFaceMeshInteractor::initInteractor for \""
                             << mesh->getName() << "\" done in " << setQTimer.getTotalTime() + solidTimer.getTotalTime()
                             << "s");
    }

    // calcForces arbeitet nicht korrekt, wenn Geo mit Bloecken
    // unterschiedlicher Leveltiefe diskretisiert -> exception
    // abfrage steht hier, weil es theoretisch sein kann, dass bei parallelen rechnungen
    // genau der block mit dem anderen level auf einem anderen prozess liegt...
    // Update: es kann u.U. passieren, dass Blöcke in der Liste nicht aktiv sin
    //(falls diese z.B. duch andere Interactoren solid gesetzt wurden)
    // diese werden nicht berücksichtigt (auch nicht beid er kraftauswertung später)
    // if( this->isRelevantForForces() )
    //{
    //   int level = -1;
    //   for( std::vector<SPtr<Block3D>>::const_iterator pos = this->transBlocks.begin(); pos!=this->transBlocks.end();
    //   ++pos)
    //      if( (*pos)->isActive() )
    //      {
    //         level = (*pos)->getLevel();
    //         break;
    //      }

    //      bool check = false;
    //      for( std::vector<SPtr<Block3D>>::const_iterator pos = this->transBlocks.begin();
    //      pos!=this->transBlocks.end(); ++pos)
    //         if( (*pos)->isActive() && (*pos)->getLevel()!=level)
    //         {
    //            (*pos)->setRank(1000);
    //            check = true;
    //            UbTupleDouble3 coords = grid.lock()->getBlockWorldCoordinates((*pos));
    //            std::cout<<(*pos)->getLevel()<<","<<(*pos)->getX1()<<","<<(*pos)->getX2()<<","<<(*pos)->getX3()<<std::endl;
    //            std::cout<<std::setprecision(15)<<val<1>(coords)<<","<<val<2>(coords)<<","<<val<3>(coords)<<std::endl<<std::endl;

    //         }
    //         if(check)
    //         {
    //            //this->grid.lock()->writeBlocks(UbStaticPathMap::getPath(UbStaticPathMap::GLOBAL)+"/error_grid",0,
    //            WbWriterVtkXmlASCII::getInstance(), false);

    //            throw UbException(UB_EXARGS,"interactor is relevant for forces,"
    //               +(std::string)" but has transblocks with different levels (wrote error_grid)"
    //               +(std::string)" -> not supportet by LBMInteractor::getForces()"
    //               +(std::string)" -> increase refineWidth");

    //         }
    //}
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::refineBlockGridToLevel(int level, double startDistance, double stopDistance)
{
    UBLOG(logDEBUG1, "D3Q27TriFaceMeshInteractor::refineBlockGridToLevel - start");

    // ToDo: evtl checken, ob man noch einen HalbraumCheck für StopDistance einbaut
    //      oder ob man schneller ist, wenn man gar keinen halbraum test macht...
    if (!grid.lock())
        throw UbException(UB_EXARGS, "Grid isn't exist!");
    if (UbMath::greater(startDistance, 0.0))
        throw UbException(UB_EXARGS, "startDistance>0.0 not supported by this interactor");
    if (UbMath::less(stopDistance, 0.0))
        throw UbException(UB_EXARGS, "stopDistance<0.0  not supported by this interactor");

    SPtr<Grid3D> bgrid    = this->grid.lock();
    GbTriFaceMesh3D &mesh = dynamic_cast<GbTriFaceMesh3D &>(*this->geoObject3D.get());

    int coarsestLevel = bgrid->getCoarsestInitializedLevel();

    std::vector<GbTriFaceMesh3D::TriFace> &triangles = *mesh.getTriangles();
    std::vector<GbTriFaceMesh3D::Vertex> &nodes      = *mesh.getNodes();

    double minX1, minX2, minX3, maxX1, maxX2, maxX3;
    float blockMinX[3], blockMaxX[3], boxCenter[3], halfBoxSize[3];
    float triPoints[3][3];

    size_t nofTriangles = (int)triangles.size();

    //#pragma omp parallel
    //#pragma omp for
    for (size_t i = 0; i < nofTriangles; i++)
    // for(int i=0; i<nofTriangles; i++)
    {
        //#pragma omp master
        //{
        //    printf_s("num_threads=%d\n", omp_get_num_threads( ));
        //}

        GbTriFaceMesh3D::TriFace &triangle = triangles[i];

        GbTriFaceMesh3D::Vertex &v1 = nodes[triangle.v1];
        GbTriFaceMesh3D::Vertex &v2 = nodes[triangle.v2];
        GbTriFaceMesh3D::Vertex &v3 = nodes[triangle.v3];

        // dreick muss normal besitzen!
        assert(!UbMath::zero(triangle.nx) || !UbMath::zero(triangle.ny) || !UbMath::zero(triangle.nz));
        // Normale muss normiert sein!
        assert((fabs(std::sqrt(triangle.nx * triangle.nx + triangle.ny * triangle.ny + triangle.nz * triangle.nz)) -
                1.0f) < 1.0E-6);

        // Halfspace um  startDistance entgegen normale verscheiebn, ansonsten werden spaeter
        // zu testende bloecke auf der dreicksrueckseite nicht getestet!!!
        GbHalfSpace3D halfSpace(
            v1.x + startDistance * triangle.nx, v1.y + startDistance * triangle.ny, v1.z + startDistance * triangle.nz,
            v2.x + startDistance * triangle.nx, v2.y + startDistance * triangle.ny, v2.z + startDistance * triangle.nz,
            v3.x + startDistance * triangle.nx, v3.y + startDistance * triangle.ny, v3.z + startDistance * triangle.nz);

        // Boundingbox um massgebliches dx erweitern -> zur Bestimmung der zu testenden Bloecke
        if (triangle.nx > 1.0E-8) {
            minX1 = triangle.getMinX(nodes) + 1.05 * triangle.nx * startDistance;
            maxX1 = triangle.getMaxX(nodes) + 1.05 * triangle.nx * stopDistance;
        } else {
            minX1 = triangle.getMinX(nodes) + 1.05 * triangle.nx * stopDistance;
            maxX1 = triangle.getMaxX(nodes) + 1.05 * triangle.nx * startDistance;
        }

        if (triangle.ny > 1.0E-8) {
            minX2 = triangle.getMinY(nodes) + 1.05 * triangle.ny * startDistance;
            maxX2 = triangle.getMaxY(nodes) + 1.05 * triangle.ny * stopDistance;
        } else {
            minX2 = triangle.getMinY(nodes) + 1.05 * triangle.ny * stopDistance;
            maxX2 = triangle.getMaxY(nodes) + 1.05 * triangle.ny * startDistance;
        }

        if (triangle.nz > 1.0E-8) {
            minX3 = triangle.getMinZ(nodes) + 1.05 * triangle.nz * startDistance;
            maxX3 = triangle.getMaxZ(nodes) + 1.05 * triangle.nz * stopDistance;
        } else {
            minX3 = triangle.getMinZ(nodes) + 1.05 * triangle.nz * stopDistance;
            maxX3 = triangle.getMaxZ(nodes) + 1.05 * triangle.nz * startDistance;
        }

        int flag = 0;
        // Levelweise alle Bloecke holen, die erweiterte BB schneiden
        // und bearbeiten
        for (int l = coarsestLevel; l < level; l++) {
            std::vector<SPtr<Block3D>> consideredBlocks;
            bgrid->getBlocksByCuboid(l, minX1, minX2, minX3, maxX1, maxX2, maxX3, consideredBlocks);
            double x1a, x2a, x3a, x1b, x2b, x3b;

            for (size_t b = 0; b < consideredBlocks.size(); b++) {
                SPtr<Block3D> block = consideredBlocks[b];
                if (block->getLevel() >= level)
                    continue;

                // start coordinaten des blocks ermitteln
                UbTupleDouble3 coords = bgrid->getBlockWorldCoordinates(block);
                UbTupleDouble3 deltas = bgrid->getBlockLengths(block);

                // Check, ob block komplett im Halbraum
                x1a = val<1>(coords);
                x1b = val<1>(coords) + val<1>(deltas);
                x2a = val<2>(coords);
                x2b = val<2>(coords) + val<2>(deltas);
                x3a = val<3>(coords);
                x3b = val<3>(coords) + val<3>(deltas);

                flag = 0;
                if (!halfSpace.ptInside(x1a, x2a, x3a))
                    flag |= (1 << 0); // 1
                if (!halfSpace.ptInside(x1b, x2a, x3a))
                    flag |= (1 << 1); // 2
                if (!halfSpace.ptInside(x1b, x2b, x3a))
                    flag |= (1 << 2); // 4
                if (!halfSpace.ptInside(x1a, x2b, x3a))
                    flag |= (1 << 3); // 8
                if (!halfSpace.ptInside(x1a, x2a, x3b))
                    flag |= (1 << 4); // 16
                if (!halfSpace.ptInside(x1b, x2a, x3b))
                    flag |= (1 << 5); // 32
                if (!halfSpace.ptInside(x1b, x2b, x3b))
                    flag |= (1 << 6); // 64
                if (!halfSpace.ptInside(x1a, x2b, x3b))
                    flag |= (1 << 7); // 128

                if (true && flag != 255) {
                    // blockseite ermitteln (skalarprodukt dreiecks-normale, vector (midTri->midCub) )
                    // je nachdem muss für den massgeblichen block start oder stopdistance verwendet werden
                    // liegt block auf pos seite -> stopdistance ansonsten startdistance
                    double skalarprod = triangle.nx * (0.5 * (x1a + x1b) - triangle.getX1Centroid(nodes)) +
                                        triangle.ny * (0.5 * (x2a + x2b) - triangle.getX2Centroid(nodes)) +
                                        triangle.nz * (0.5 * (x3a + x3b) - triangle.getX3Centroid(nodes));

                    double blockdelta = 1.05 * stopDistance;
                    if (skalarprod < 1.E-8)
                        blockdelta = -1.05 * startDistance; // startDistance<0!!
                    else if (fabs(skalarprod) < 1.E-8)
                        blockdelta = 1.05 * UbMath::max(-startDistance, stopDistance);

                    // block anpassen
                    blockMinX[0] = (float)(val<1>(coords) - blockdelta);
                    blockMinX[1] = (float)(val<2>(coords) - blockdelta);
                    blockMinX[2] = (float)(val<3>(coords) - blockdelta);

                    blockMaxX[0] = (float)(val<1>(coords) + val<1>(deltas) + blockdelta);
                    blockMaxX[1] = (float)(val<2>(coords) + val<2>(deltas) + blockdelta);
                    blockMaxX[2] = (float)(val<3>(coords) + val<3>(deltas) + blockdelta);

                    boxCenter[0] = (float)(0.5 * (blockMaxX[0] + blockMinX[0]));
                    boxCenter[1] = (float)(0.5 * (blockMaxX[1] + blockMinX[1]));
                    boxCenter[2] = (float)(0.5 * (blockMaxX[2] + blockMinX[2]));

                    halfBoxSize[0] = (float)(0.5 * (blockMaxX[0] - blockMinX[0]));
                    halfBoxSize[1] = (float)(0.5 * (blockMaxX[1] - blockMinX[1]));
                    halfBoxSize[2] = (float)(0.5 * (blockMaxX[2] - blockMinX[2]));

                    GbTriFaceMesh3D::Vertex &v1_ = nodes[triangle.v1];
                    GbTriFaceMesh3D::Vertex &v2_ = nodes[triangle.v2];
                    GbTriFaceMesh3D::Vertex &v3_ = nodes[triangle.v3];

                    triPoints[0][0] = v1_.x;
                    triPoints[0][1] = v1_.y;
                    triPoints[0][2] = v1_.z;
                    triPoints[1][0] = v2_.x;
                    triPoints[1][1] = v2_.y;
                    triPoints[1][2] = v2_.z;
                    triPoints[2][0] = v3_.x;
                    triPoints[2][1] = v3_.y;
                    triPoints[2][2] = v3_.z;

                    // wenn block dreick schneidet, dann muss er verfeinert werden
                    if (GbMeshTools3D::triBoxOverlap(boxCenter, halfBoxSize, triPoints)) {
                        bgrid->expandBlock(block->getX1(), block->getX2(), block->getX3(), block->getLevel());
                    }
                }
            }
        }
    }
    UBLOG(logDEBUG1, " - refine done");
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::updateMovedGeometry(const real &timeStep) {}
////////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::recursiveGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys,
                                                   const short &zs, const FLAGS &type)
{
    // Algorithmus zum Füllen eines Polyeders, ausgehend vom Saatpunkt xs,ys,zs

    // Saatknoten einfärben
    if (flagfield(xs, ys, zs) == UNDEF_FLAG) {
        flagfield(xs, ys, zs) = type;

        if (flagfield.indicesInRange(xs + 1, ys, zs))
            this->recursiveGridFill(flagfield, xs + 1, ys, zs, type);
        if (flagfield.indicesInRange(xs, ys + 1, zs))
            this->recursiveGridFill(flagfield, xs, ys + 1, zs, type);
        if (flagfield.indicesInRange(xs, ys, zs + 1))
            this->recursiveGridFill(flagfield, xs, ys, zs + 1, type);
        if (flagfield.indicesInRange(xs - 1, ys, zs))
            this->recursiveGridFill(flagfield, xs - 1, ys, zs, type);
        if (flagfield.indicesInRange(xs, ys - 1, zs))
            this->recursiveGridFill(flagfield, xs, ys - 1, zs, type);
        if (flagfield.indicesInRange(xs, ys, zs - 1))
            this->recursiveGridFill(flagfield, xs, ys, zs - 1, type);
    }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::iterativeGridFill(CbArray3D<FLAGS> &flagfield, const short &xs, const short &ys,
                                                   const short &zs, const FLAGS &type)
{
    std::stack<UbTupleInt3> stck;
    stck.push(UbTupleInt3(xs, ys, zs));

    int x, y, z;

    while (!stck.empty()) {
        x = val<1>(stck.top());
        y = val<2>(stck.top());
        z = val<3>(stck.top());
        stck.pop();

        FLAGS &flagType = flagfield(x, y, z);

        if (flagType == UNDEF_FLAG) {
            flagType = type;

            if (flagfield.indicesInRange(x + 1, y, z))
                stck.push(UbTupleInt3(x + 1, y, z));
            if (flagfield.indicesInRange(x, y + 1, z))
                stck.push(UbTupleInt3(x, y + 1, z));
            if (flagfield.indicesInRange(x, y, z + 1))
                stck.push(UbTupleInt3(x, y, z + 1));
            if (flagfield.indicesInRange(x - 1, y, z))
                stck.push(UbTupleInt3(x - 1, y, z));
            if (flagfield.indicesInRange(x, y - 1, z))
                stck.push(UbTupleInt3(x, y - 1, z));
            if (flagfield.indicesInRange(x, y, z - 1))
                stck.push(UbTupleInt3(x, y, z - 1));
        }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 D3Q27TriFaceMeshInteractor::getForces()
{
    // FeTriFaceMesh3D* feMesh = dynamic_cast<FeTriFaceMesh3D*>(this->geoObject3D.get());
    // if(!feMesh)
    //{
    //   return D3Q19AMRInteractor::getForces();
    //}
    ////return getForcesTriangle();
    // this->calculateForces();

    real forceX1 = 0.0;
    real forceX2 = 0.0;
    real forceX3 = 0.0;

    // double area = 0.0;

    // vector<FeTriFaceMesh3D::VertexAttributes>* attributes = feMesh->getAttributes();

    // for(size_t i=0; i<attributes->size(); i++)
    //{
    //   FeTriFaceMesh3D::VertexAttributes& attribut = (*attributes)[i];
    //   area = attribut.getArea();
    //   forceX1 += attribut.getFX()*area;
    //   forceX2 += attribut.getFY()*area;
    //   forceX3 += attribut.getFZ()*area;
    //}
    return { forceX1, forceX2, forceX3 };
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 D3Q27TriFaceMeshInteractor::getForcesTriangle()
{
    real forceX1 = 0.0;
    real forceX2 = 0.0;
    real forceX3 = 0.0;

    // D3Q19BlockGrid& grid.lock() = dynamic_cast<D3Q19BlockGrid&>(*this->grid.lock());
    ////   CoordinateTransformation3D *trafo = this->grid.lock()->getTransformation();
    ////   int minLevel = this->grid.lock()->getFinestInitializedLevel();
    ////    double scaleX = trafo->getX1CoordinateScaling()/(1<<minLevel);
    ////    double scaleY = trafo->getX2CoordinateScaling()/(1<<minLevel);
    ////    double scaleZ = trafo->getX3CoordinateScaling()/(1<<minLevel);
    ////    int blocknx1 = grid.lock()->getBlockNX1();
    ////    int blocknx2 = grid.lock()->getBlockNX2();
    ////    int blocknx3 = grid.lock()->getBlockNX3();
    ////    double xOffset = trafo->getX1CoordinateOffset();
    ////    double yOffset = trafo->getX2CoordinateOffset();
    ////    double zOffset = trafo->getX3CoordinateOffset();
    // vector<D3Q19Real> collFactors = ((D3Q19Calculator*)grid.lock()->getCalculator())->getCollisionsFactors();

    ////for (int i=0;i<(int)gbTriangle3DInteractors.size(); i++)
    ////{
    ////   GbTriangle3D* tri = (GbTriangle3D*)gbTriangle3DInteractors[i]->getGbObject3D();

    ////   double px0 = tri->getX1Centroid();
    ////   double py0 = tri->getX2Centroid();
    ////   double pz0 = tri->getX3Centroid();
    ////   double px = px0-xOffset;
    ////   double py = py0-yOffset;
    ////   double pz = pz0-zOffset;
    ////   px = px/scaleX;
    ////   py = py/scaleY;
    ////   pz = pz/scaleZ;
    ////   int x1 = (int)px;
    ////   int y1 = (int)py;
    ////   int z1 = (int)pz;
    ////   AMR3DBlock* block = this->grid.lock()->getBlock(x1,y1,z1,minLevel);
    ////   if(!block)  block = this->grid.lock()->getSuperBlock(x1,y1,z1,minLevel);
    ////   if(!block) throw UbException(__FILE__,__LINE__,"kein Block ...");

    ////   double collFactor = collFactors[block->getLevel()];
    ////   double nodeDistance = grid.lock()->getNodeDeltaX(block->getLevel());
    ////   double bertX1 = ((px0-xOffset)/(1000.*nodeDistance));
    ////   double bertX2 = ((py0-yOffset)/(1000.*nodeDistance));
    ////   double bertX3 = ((pz0-zOffset)/(1000.*nodeDistance));
    ////   int abstaendeX1 = (int)(bertX1*1000.);
    ////   int abstaendeX2 = (int)(bertX2*1000.);
    ////   int abstaendeX3 = (int)(bertX3*1000.);
    ////   int posW = abstaendeX1 - block->getX1Index()*blocknx1;
    ////   int posS = abstaendeX2 - block->getX2Index()*blocknx2;
    ////   int posB = abstaendeX3 - block->getX3Index()*blocknx3;
    ////   int posE=posW+1;
    ////   int posN=posS+1;
    ////   int posT=posB+1;

    ////   D3Q19BlockDescriptor *bvd = dynamic_cast<D3Q19BlockDescriptor*>(block->getBlockDescriptor());
    ////   if(!bvd) throw UbException(__FILE__,__LINE__,"kein Bvd ...");

    ////   CbUniformMatrix4D<double,IndexerX1X2X3X4>* tempdistributions = bvd->getTempDistributionMatrix();
    ////   D3Q19BCMatrix<D3Q19BoundaryCondition> *bcMatrix = bvd->getBcMatrix();

    ////   UbTupleDouble6 stresses;
    ////   double dX = px0-this->geoObject3D->getX1Centroid();
    ////   double dY = py0-this->geoObject3D->getX2Centroid();
    ////   double dZ = pz0-this->geoObject3D->getX3Centroid();
    ////   if(dX<=0.0 && dY<=0.0 && dZ<=0.0)
    ////   {
    ////      double *fWSB  = tempdistributions->getStartAdressOfSortedArray(posW,posS,posB,0);
    ////      if(bcMatrix->isFluid(posW,posS,posB)) stresses = D3Q19System::getIncompStresses(fWSB, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX<=0.0 && dY>0.0 && dZ<=0.0)
    ////   {
    ////      double *fWNB  = tempdistributions->getStartAdressOfSortedArray(posW,posN,posB,0);
    ////      if(bcMatrix->isFluid(posW,posN,posB)) stresses = D3Q19System::getIncompStresses(fWNB, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX<=0.0 && dY<=0.0 && dZ>0.0)
    ////   {
    ////      double *fWST  = tempdistributions->getStartAdressOfSortedArray(posW,posS,posT,0);
    ////      if(bcMatrix->isFluid(posW,posS,posT)) stresses = D3Q19System::getIncompStresses(fWST, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX<=0.0 && dY>0.0 && dZ>0.0)
    ////   {
    ////      double *fWNT  = tempdistributions->getStartAdressOfSortedArray(posW,posN,posT,0);
    ////      if(bcMatrix->isFluid(posW,posN,posT)) stresses = D3Q19System::getIncompStresses(fWNT, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX>0.0 && dY<=0.0 && dZ<=0.0)
    ////   {
    ////      double *fESB  = tempdistributions->getStartAdressOfSortedArray(posE,posS,posB,0);
    ////      if(bcMatrix->isFluid(posE,posS,posB)) stresses = D3Q19System::getIncompStresses(fESB, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX>0.0 && dY>0.0 && dZ<=0.0)
    ////   {
    ////      double *fENB  = tempdistributions->getStartAdressOfSortedArray(posE,posN,posB,0);
    ////      if(bcMatrix->isFluid(posE,posN,posB)) stresses = D3Q19System::getIncompStresses(fENB, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX>0.0 && dY<=0.0 && dZ>0.0)
    ////   {
    ////      double *fEST  = tempdistributions->getStartAdressOfSortedArray(posE,posS,posT,0);
    ////      if(bcMatrix->isFluid(posE,posS,posT)) stresses = D3Q19System::getIncompStresses(fEST, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else if(dX>0.0 && dY>0.0 && dZ>0.0)
    ////   {
    ////      double *fENT  = tempdistributions->getStartAdressOfSortedArray(posE,posN,posT,0);
    ////      if(bcMatrix->isFluid(posE,posN,posT)) stresses = D3Q19System::getIncompStresses(fENT, collFactor );
    ////      else cout<<__LINE__<<" nicht fluid ...";
    ////   }
    ////   else cout<<"punkt mit:"<<dX<<" "<<dY<<" "<<dZ<<" ist nicht bei \n";

    ////   double S11 = val<1>(stresses);
    ////   double S22 = val<2>(stresses);
    ////   double S33 = val<3>(stresses);
    ////   double S12 = val<4>(stresses);
    ////   double S13 = val<5>(stresses);
    ////   double S23 = val<6>(stresses);

    ////   GbVector3D normal = tri->getNormal();
    ////   double nx = normal.X1();
    ////   double ny = normal.X2();
    ////   double nz = normal.X3();
    ////   double area = tri->getArea();

    ////   double Fx1 = area*(S11*nx+S12*ny+S13*nz);
    ////   double Fy1 = area*(S12*nx+S22*ny+S23*nz);
    ////   double Fz1 = area*(S13*nx+S23*ny+S33*nz);
    ////   forceX1 += Fx1;
    ////   forceX2 += Fy1;
    ////   forceX3 += Fz1;
    ////}
    return { forceX1, forceX2, forceX3 };
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::calculateForces()
{
    // FeTriFaceMesh3D* feMesh = dynamic_cast<FeTriFaceMesh3D*>(this->geoObject3D.get());
    // if(!feMesh) throw UbException(UB_EXARGS,"geoObject is not a FeTriFaceMesh3D!");

    // if(this->stressMode == STRESSNORMAL) this->calculateStresses();
    // else if(this->stressMode == STRESSALTERNATIV) this->calculateStressesAlternativ();

    // vector<FeTriFaceMesh3D::VertexAttributes>* attributes = feMesh->getAttributes();

    // for (int i=0;i<(int)attributes->size() ;i++)
    //{
    //   FeTriFaceMesh3D::VertexAttributes& attribut = (*attributes)[i];
    //   attribut.setFX(0.0);
    //   attribut.setFY(0.0);
    //   attribut.setFZ(0.0);
    //   attribut.setArea(0.0);
    //}
    // vector<GbTriFaceMesh3D::TriFace>& triangles = *feMesh->getTriangles();
    // vector<GbTriFaceMesh3D::Vertex>&  nodes = *feMesh->getNodes();
    // for (size_t i=0; i<triangles.size(); i++)
    //{
    //   GbTriFaceMesh3D::TriFace& triangle = triangles[i];
    //   FeTriFaceMesh3D::VertexAttributes& vAttribut1 = (*attributes)[triangle.v1];
    //   FeTriFaceMesh3D::VertexAttributes& vAttribut2 = (*attributes)[triangle.v2];
    //   FeTriFaceMesh3D::VertexAttributes& vAttribut3 = (*attributes)[triangle.v3];
    //   UbTupleDouble6& stressesP1 = vAttribut1.getStresses();
    //   UbTupleDouble6& stressesP2 = vAttribut2.getStresses();
    //   UbTupleDouble6& stressesP3 = vAttribut3.getStresses();
    //   double p1S11 = val<1>(stressesP1); double p2S11 = val<1>(stressesP2); double p3S11 = val<1>(stressesP3);
    //   double p1S22 = val<2>(stressesP1); double p2S22 = val<2>(stressesP2); double p3S22 = val<2>(stressesP3);
    //   double p1S33 = val<3>(stressesP1); double p2S33 = val<3>(stressesP2); double p3S33 = val<3>(stressesP3);
    //   double p1S12 = val<4>(stressesP1); double p2S12 = val<4>(stressesP2); double p3S12 = val<4>(stressesP3);
    //   double p1S13 = val<5>(stressesP1); double p2S13 = val<5>(stressesP2); double p3S13 = val<5>(stressesP3);
    //   double p1S23 = val<6>(stressesP1); double p2S23 = val<6>(stressesP2); double p3S23 = val<6>(stressesP3);

    //   triangle.calculateNormal(nodes);
    //   double nx = triangle.nx;
    //   double ny = triangle.ny;
    //   double nz = triangle.nz;
    //   double area = 0.3333*triangle.getArea(nodes);

    //   if(UbMath::lessEqual(area,0.0)) cout<<__FILE__<<" "<<__LINE__<<" area <= 0 "<<endl;

    //   double Fx1 =
    //   area*(0.333*(p1S11*nx+p1S12*ny+p1S13*nz)+0.333*(p2S11*nx+p2S12*ny+p2S13*nz)+0.333*(p3S11*nx+p3S12*ny+p3S13*nz));
    //   double Fx2 = Fx1;
    //   double Fx3 = Fx1;

    //   double Fy1 =
    //   area*(0.333*(p1S12*nx+p1S22*ny+p1S23*nz)+0.333*(p2S12*nx+p2S22*ny+p2S23*nz)+0.333*(p3S12*nx+p3S22*ny+p3S23*nz));
    //   double Fy2 = Fy1;
    //   double Fy3 = Fy1;

    //   double Fz1 =
    //   area*(0.333*(p1S13*nx+p1S23*ny+p1S33*nz)+0.333*(p2S13*nx+p2S23*ny+p2S33*nz)+0.333*(p3S13*nx+p3S23*ny+p3S33*nz));
    //   double Fz2 = Fz1;
    //   double Fz3 = Fz1;
    //   //  cout<<Fx1<<" "<<Fy1<<" "<<Fz1<<endl;
    //   vAttribut1.addFX(Fx1);    vAttribut2.addFX(Fx2);    vAttribut3.addFX(Fx3);
    //   vAttribut1.addFY(Fy1);    vAttribut2.addFY(Fy2);    vAttribut3.addFY(Fy3);
    //   vAttribut1.addFZ(Fz1);    vAttribut2.addFZ(Fz2);    vAttribut3.addFZ(Fz3);
    //   vAttribut1.addArea(area); vAttribut2.addArea(area); vAttribut3.addArea(area);
    //}
    // for (size_t i=0; i<attributes->size(); i++)
    //{
    //   FeTriFaceMesh3D::VertexAttributes& attribut = (*attributes)[i];

    //   double newFX = attribut.getFX()/attribut.getArea();
    //   double newFY = attribut.getFY()/attribut.getArea();
    //   double newFZ = attribut.getFZ()/attribut.getArea();
    //   //if(i==100) cout<<"F:"<<newFX<<" "<<newFY<<" "<<newFZ<<endl;
    //   //double oldFX = p->getOldFX();
    //   //double oldFY = p->getOldFY();
    //   //int alphaSteps = p->getFilteringSteps();
    //   //double alpha = 1.0;
    //   //if(alphaSteps != 0)
    //   //{
    //   //   alpha = (1.0-alphaSteps*0.1);
    //   // //  cout<<p->toString()<<" alpha:"<<alpha<<" steps:"<<alphaSteps<<endl;
    //   //   p->reduceFilteringSteps();
    //   //}
    //   //newFX = (1.-alpha)*oldFX+alpha*newFX;
    //   //newFY = (1.-alpha)*oldFY+alpha*newFY;

    //   attribut.setFX(newFX);
    //   attribut.setFY(newFY);
    //   attribut.setFZ(newFZ);
    //   //cout<<i<<" "<<newFX<<" "<<newFY<<" "<<newFZ<<endl;
    //   //cout<<i<<" "<<p->toString()<<endl;

    //}
}
//////////////////////////////////////////////////////////////////////////
string D3Q27TriFaceMeshInteractor::toString()
{
    stringstream ss;
    ss << "D3Q27TriFaceMeshInteractor[label=D3Q27TriFaceMeshInteractor";
    if (this->isSolid())
        ss << ", solid";
    if (this->isInverseSolid())
        ss << ", inversesolid";
    if (this->isTimeDependent())
        ss << ", timedependent";
    if (geoObject3D != NULL)
        ss << ", AMR3DInteractor: " << geoObject3D->toString();
    ss << "]";

    return ss.str();
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::reinitWithStoredQs(const real & /*timeStep*/)
{
    // alle solid Bloecke wieder solid setzen
    std::vector<SPtr<Block3D>> &solidBlocks = this->getSolidBlockSet();
    for (size_t i = 0; i < solidBlocks.size(); i++) {
        solidBlocks[i]->setActive(false); //<- quick n dirty
    }

    // alle solid-nodes wieder solid setzen (solids die quasi in den TransBloecken liegen)
    std::map<SPtr<Block3D>, std::set<UbTupleInt3>>::iterator it1;
    for (it1 = this->solidNodeIndicesMap.begin(); it1 != this->solidNodeIndicesMap.end(); ++it1) {
        SPtr<Block3D> block = it1->first;

        SPtr<ILBMKernel> kernel           = block->getKernel();
        SPtr<BCArray3D> bcMatrix          = kernel->getBCSet()->getBCArray();
        std::set<UbTupleInt3> &indicesSet = it1->second;

        for (std::set<UbTupleInt3>::iterator setIt = indicesSet.begin(); setIt != indicesSet.end(); ++setIt) {
            bcMatrix->setSolid(val<1>(*setIt), val<2>(*setIt), val<3>(*setIt));
        }
    }

    // BCS WIEDERHERSTELLEN
    std::map<SPtr<Block3D>, std::map<UbTupleInt3, std::vector<float>>>::iterator it;
    for (it = bcNodeIndicesAndQsMap.begin(); it != bcNodeIndicesAndQsMap.end(); ++it) {
        SPtr<Block3D> block      = it->first;
        SPtr<ILBMKernel> kernel  = block->getKernel();
        SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

        std::map<UbTupleInt3, std::vector<float>>::iterator it2;
        for (it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            const UbTupleInt3 &pos = it2->first;
            std::vector<float> qs  = it2->second;

            // SG_27.08.2010
            if (bcMatrix->isSolid(val<1>(pos), val<2>(pos), val<3>(pos)))
                continue;

            SPtr<BoundaryConditions> bc = bcMatrix->getBC(val<1>(pos), val<2>(pos), val<3>(pos));
            if (!bc) {
                bc = SPtr<BoundaryConditions>(new BoundaryConditions);
                bcMatrix->setBC(val<1>(pos), val<2>(pos), val<3>(pos), bc);
            }

            double x1w = qs[D3Q27System::FENDDIR + 1 + 0];
            double x2w = qs[D3Q27System::FENDDIR + 1 + 1];
            double x3w = qs[D3Q27System::FENDDIR + 1 + 2];

            // TODO: HACK GEHOERT NICHT HIERHIER!!! - start
            // es handelt sich un ein statisches Objekt und beim Propeller gibt
            // es Schwierigkeiten an den Flügelspitzen, dass kann daher kommen,
            // dass dort zuviel bc-flaggs sind und mit Geschwindigkeit ergibt dies ziemlich grosse Werte
            bc->setBoundaryVelocityX1(0.0);
            bc->setBoundaryVelocityX2(0.0);
            bc->setBoundaryVelocityX3(0.0);
            // TODO: HACK GEHOERT NICHT HIERHIER!!! - end

            bool gotQs = false;
            for (int fdir = D3Q27System::FSTARTDIR; fdir <= D3Q27System::FENDDIR; fdir++) {
                if (UbMath::greater(qs[fdir], -1.0) && UbMath::less(qs[fdir], bc->getQ(fdir))) {
                    gotQs = true;
                    for (size_t index = 0; index < this->BCs.size(); index++)
                        this->BCs[index]->adaptBCForDirection(*this, bc, x1w, x2w, x3w, qs[fdir], fdir);
                }
            }

            if (gotQs)
                for (size_t index = 0; index < this->BCs.size(); index++)
                    this->BCs[index]->adaptBC(*this, bc, x1w, x2w, x3w);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void D3Q27TriFaceMeshInteractor::updateInteractor(const real &timestep)
{
    D3Q27Interactor::updateInteractor(timestep);
}
