//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup geometry3d
//! \ingroup basics
//! \{
//! \author Konstantin Kutscher, Soeren Textor, Sebastian Geller
//=======================================================================================
#ifndef GBVOXELMATRIX3D_H
#define GBVOXELMATRIX3D_H

#include <cmath>
#include <vector>

#include <basics/container/CbArray3D.h>
#include <basics/utilities/UbObserver.h>
#include <geometry3d/GbObject3D.h>

#include <PointerDefinitions.h>

class GbLine3D;
class GbTriangle3D;
class GbObject3DCreator;

class GbVoxelMatrix3D : public GbObject3D, public UbObserver
{
public:
    using Matrix3D = CbArray3D<float>;
    static const float SOLID;
    static const float FLUID;
    enum Endian { BigEndian, LittleEndian };

    GbVoxelMatrix3D();
    GbVoxelMatrix3D(int nx1, int nx2, int nx3, float initVal, double lowerThreshold = 0, double upperThreshold = 0);
    GbVoxelMatrix3D(const Matrix3D &voxelMatrix, double lowerThreshold = 0, double upperThreshold = 0);
    ~GbVoxelMatrix3D() override = default;

    void finalize() override{};
    GbVoxelMatrix3D *clone() override;

    /*=======================================================================*/
    Matrix3D::reference operator()(const Matrix3D::size_type &x1, const Matrix3D::size_type &x2,
                                   const Matrix3D::size_type &x3)
    {
        return voxelMatrix(x1, x2, x3);
    }
    /*=======================================================================*/
    Matrix3D::const_reference operator()(const Matrix3D::size_type &x1, const Matrix3D::size_type &x2,
                                         const Matrix3D::size_type &x3) const
    {
        return voxelMatrix(x1, x2, x3);
    }
    /*=======================================================================*/
    void setTransferViaFilename(bool transferViaFilename, std::string filename)
    {
        this->filename            = filename;
        this->transferViaFilename = transferViaFilename;
    }
    void setThreshold(double lowerThreshold, double upperThreshold)
    {
        this->lowerThreshold = lowerThreshold;
        this->upperThreshold = upperThreshold;
    }
    void setAddSurfaceTriangleSetFlag(bool flag) { this->addSurfaceTriangleSetFlag = flag; }

    /*=======================================================================*/
    void setVoxelMatrixMininum(double minX1, double minX2, double minX3)
    {
        this->minX1 = minX1;
        this->minX2 = minX2;
        this->minX3 = minX3;
    }
    void setVoxelMatrixMinX1(double minX1) { this->minX1 = minX1; }
    void setVoxelMatrixMinX2(double minX2) { this->minX2 = minX2; }
    void setVoxelMatrixMinX3(double minX3) { this->minX3 = minX3; }

    /*=======================================================================*/
    void setVoxelMatrixDelta(double deltaX1, double deltaX2, double deltaX3)
    {
        this->deltaX1 = deltaX1;
        this->deltaX2 = deltaX2;
        this->deltaX3 = deltaX3;
    }
    void setVoxelMatrixDeltaX1(double deltaX1) { this->deltaX1 = deltaX1; }
    void setVoxelMatrixDeltaX2(double deltaX2) { this->deltaX2 = deltaX2; }
    void setVoxelMatrixDeltaX3(double deltaX3) { this->deltaX3 = deltaX3; }

    /*=======================================================================*/
    double getX1Centroid() override { return 0.5 * (minX1 + this->getX1Maximum()); }
    double getX1Minimum() override { return minX1; }
    double getX1Maximum() override { return minX1 + deltaX1 * voxelMatrix.getNX1(); }

    double getX2Centroid() override { return 0.5 * (minX2 + this->getX2Maximum()); }
    double getX2Minimum() override { return minX2; }
    double getX2Maximum() override { return minX2 + deltaX2 * voxelMatrix.getNX2(); }

    double getX3Centroid() override { return 0.5 * (this->getX3Minimum() + this->getX3Maximum()); }
    double getX3Minimum() override { return minX3; }
    double getX3Maximum() override { return minX3 + deltaX3 * voxelMatrix.getNX3(); }

    double getLengthX1() { return this->getX1Maximum() - minX1; }
    double getLengthX2() { return this->getX2Maximum() - minX2; }
    double getLengthX3() { return this->getX3Maximum() - minX3; }

    void setCenterCoordinates(const double &x1, const double &x2, const double &x3) override;
    void translate(const double &tx1, const double &tx2, const double &tx3) override;

    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p, bool &pointIsOnBoundary) override;
    bool isPointInGbObject3D(const double &x1p, const double &x2p, const double &x3p) override;
    bool isCellInsideGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                const double &x2b, const double &x3b) override;
    bool isCellCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                 const double &x2b, const double &x3b) override;
    bool isCellInsideOrCuttingGbObject3D(const double &x1a, const double &x2a, const double &x3a, const double &x1b,
                                         const double &x2b, const double &x3b) override;
    // double getCellVolumeInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double&
    // x1b,const double& x2b,const double& x3b);

    GbPoint3D *calculateInterSectionPoint3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/)
    {
        throw UbException(__FILE__, __LINE__, UB_FUNCTION, "not implemented");
    }
    GbLine3D *createClippedLine3D(GbPoint3D & /*point1*/, GbPoint3D & /*point2*/) override
    {
        throw UbException(__FILE__, __LINE__, UB_FUNCTION, "not implemented");
    }

    std::vector<GbTriangle3D *> getSurfaceTriangleSet() override;
    void addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles) override;

    bool hasRaytracing() override { return true; }
    /*|r| must be 1! einheitsvector!!*/
    double getIntersectionRaytraceFactor(const double &x1, const double &x2, const double &x3, const double &rx1,
                                         const double &rx2, const double &rx3) override;

    std::string toString() override;

    // virtuelle Methoden von UbObserver
    void objectChanged(UbObservable *changedObject) override {}
    void objectWillBeDeleted(UbObservable *objectForDeletion) override {}

    template <class T>
    void readMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian);
    template <class T>
    void readBufferedMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian);
    void readMatrixFromVtiASCIIFile(std::string filename);

    void rotate90aroundX();
    void rotate90aroundY();
    void rotate90aroundZ();
    void rotate90aroundX(double cX1, double cX2, double cX3);
    void rotate90aroundY(double cX1, double cX2, double cX3);
    void rotate90aroundZ(double cX1, double cX2, double cX3);
    void mirrorX();
    void mirrorY();
    void mirrorZ();

    void rotateAroundY(double theta);

    void writeToLegacyVTKASCII(const std::string &fileName);
    void writeToLegacyVTKBinary(const std::string &fileName);
    void writeToVTKImageDataASCII(const std::string &fileName);
    void writeToVTKImageDataAppended(const std::string &fileName);

    void setClosedVoidSpaceToSolid();

    void calculateNumberOfSolidAndFluid();
    long getNumberOfSolid();
    long getNumberOfFluid();

protected:
    void findFluidNeighbor(int cx1, int cx2, int cx3);
    Matrix3D voxelMatrixTemp;
    CbArray3D<char> flagMatrix;

    std::vector<int> x1Nbr;
    std::vector<int> x2Nbr;
    std::vector<int> x3Nbr;

    std::vector<int> x1NbrTemp;
    std::vector<int> x2NbrTemp;
    std::vector<int> x3NbrTemp;

    using GbObject3D::isPointInGbObject3D; // Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht
                                           // ausprogrammieren, welche sonst hier "ueberdeckt" waere

protected:
    // for transfer
    std::string filename;
    bool transferViaFilename{ false };

    bool addSurfaceTriangleSetFlag{ true };

    int nodesX1{ 0 };
    int nodesX2{ 0 };
    int nodesX3{ 0 };
    double lowerThreshold{ 0.0 }, upperThreshold{ 0.0 };

    double minX1{ 0.0 };
    double minX2{ 0.0 };
    double minX3{ 0.0 };
    double deltaX1{ 1.0 };
    double deltaX2{ 1.0 };
    double deltaX3{ 1.0 };

    Matrix3D voxelMatrix;

    long numberOfSolid;
    long numberOfFluid;
};

//////////////////////////////////////////////////////////////////////////
template <class T>
void GbVoxelMatrix3D::readMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian)
{
    using namespace std;
    // UBLOG(logINFO,"GbVoxelMatrix3D::readMatrixFromFile \""<<filename<<"\"
    // nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - start");
    ifstream in(filename.c_str(), ios::binary);
    if (!in)
        throw UbException(UB_EXARGS, "could not open file " + filename);

    in.seekg(0, ios::end);                 // Ende springen
    fstream::off_type length = in.tellg(); // Position abfragen
    in.seekg(0, ios::beg);                 // An den Anfang springen

    // UBLOG(logINFO,"number of nodes = "<<nodesX1*nodesX2*nodesX3*sizeof(T)<<" file size = "<<(long)length);
    // if( (nodesX1*nodesX2*nodesX3)*sizeof(float) != (long)length )
    unsigned long long nofn = (unsigned long long)nodesX1 * (unsigned long long)nodesX2 * (unsigned long long)nodesX3 *
                              (unsigned long long)sizeof(T);
    if (nofn != (unsigned long long)length) {
        throw UbException(UB_EXARGS, "number of nodes(" + UbSystem::toString(nofn) + ") doesn't match file size(" +
                                         UbSystem::toString((long)length) + ")");
    }

    // UBLOG(logINFO,"  - create GbVoxelMatrix3D");
    // GbVoxelMatrix3D* voxelGeo = new GbVoxelMatrix3D(nodesX1,nodesX2,nodesX3,GbVoxelMatrix3D::FLUID);
    voxelMatrix = Matrix3D(nodesX1, nodesX2, nodesX3, GbVoxelMatrix3D::FLUID);

    // UBLOG(logINFO,"  - init values");
    // float val;
    T val;
    for (int x3 = 0; x3 < nodesX3; x3++)
        for (int x2 = 0; x2 < nodesX2; x2++)
            for (int x1 = 0; x1 < nodesX1; x1++) {
                // in.read((char*)&val,sizeof(float));
                in.read((char *)&val, sizeof(T));
                if (endian == BigEndian)
                    UbSystem::swapByteOrder((unsigned char *)(&(val)), sizeof(T));
                // if( UbMath::equal((double)val, threshold) )
                // if( UbMath::greater((double)val, threshold) )
                if ((double)val >= lowerThreshold && (double)val <= upperThreshold) {
                    (voxelMatrix)(x1, x2, x3) = GbVoxelMatrix3D::SOLID;
                }
                //(voxelMatrix)(x1, x2, x3) = (float)val;
            }

    // UBLOG(logINFO,"GbVoxelMatrix3D::readMatrixFromFile \""<<filename<<"\"
    // nodes("<<nodesX1<<"/"<<nodesX2<<"/"<<nodesX3<<") - end");
}

//////////////////////////////////////////////////////////////////////////
template <class T>
void GbVoxelMatrix3D::readBufferedMatrixFromRawFile(std::string filename, GbVoxelMatrix3D::Endian endian)
{
    using namespace std;
    UBLOG(logINFO, "GbVoxelMatrix3D::readMatrixFromRawFile \"" << filename << "\" nodes(" << nodesX1 << "/" << nodesX2
                                                               << "/" << nodesX3 << ") - start");

    FILE *file;
    file = fopen(filename.c_str(), "rb");
    if (file == NULL) {
        throw UbException(UB_EXARGS, "Could not open file " + filename);
    }

    // obtain file size:
    fseek(file, 0, SEEK_END);
    unsigned long int length = ftell(file);
    rewind(file);

    UBLOG(logINFO, "number of nodes = " << (long)nodesX1 * (long)nodesX2 * (long)nodesX3 << " file size = " << length);

    unsigned long int nofn = (long)nodesX1 * (long)nodesX2 * (long)nodesX3 * (long)sizeof(T);
    if (nofn != length) {
        // throw UbException(UB_EXARGS, "number of nodes("+UbSystem::toString(nofn)+") doesn't match file
        // size("+UbSystem::toString(length)+")");
    }

    UBLOG(logINFO, "  - create GbVoxelMatrix3D");
    voxelMatrix = Matrix3D(nodesX1, nodesX2, nodesX3, GbVoxelMatrix3D::FLUID);

    CbArray3D<T> readMatrix(nodesX1, nodesX2, nodesX3);

    UBLOG(logINFO, "  - read file to matrix");
    fread(readMatrix.getStartAdressOfSortedArray(0, 0, 0), sizeof(T), readMatrix.getDataVector().size(), file);
    fclose(file);

    UBLOG(logINFO, "  - init values");

    numberOfSolid = 0;
    T val;
    for (int x3 = 0; x3 < nodesX3; x3++)
        for (int x2 = 0; x2 < nodesX2; x2++)
            for (int x1 = 0; x1 < nodesX1; x1++) {
                val = readMatrix(x1, x2, x3);

                if (endian == BigEndian) {
                    UbSystem::swapByteOrder((unsigned char *)(&(val)), sizeof(T));
                }

                if ((double)val >= lowerThreshold && (double)val <= upperThreshold) {
                    voxelMatrix(x1, x2, x3) = GbVoxelMatrix3D::SOLID;
                }
            }

    UBLOG(logINFO, "GbVoxelMatrix3D::readMatrixFromRawFile \"" << filename << "\" nodes(" << nodesX1 << "/" << nodesX2
                                                               << "/" << nodesX3 << ") - end");
}

//#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
// UB_AUTO_RUN_NAMED(SF::registerType<GbVoxelMatrix3D>("GbVoxelMatrix3D"), SF_GbVoxelMatrix3D);
// UB_AUTO_RUN_NAMED((SF::registerBaseAndDerived< GbObject3D, GbVoxelMatrix3D >()), SF_GbVoxelMatrix3D_BD1);
// UB_AUTO_RUN_NAMED((SF::registerBaseAndDerived< UbObserver, GbVoxelMatrix3D>()), SF_GbVoxelMatrix3D_BD2);
//#endif //RCF_USE_SF_SERIALIZATION

#endif

//! \}
