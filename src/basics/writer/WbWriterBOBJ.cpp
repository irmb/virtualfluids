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
//! \addtogroup writer
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifdef CAB_ZLIB
#include <basics/utilities/UbLogger.h>
#include <basics/writer/WbWriterBOBJ.h>
#include <cstring>

#include <zlib.h>

using namespace std;
/*===============================================================================*/
std::string WbWriterBOBJ::writeTriangles(const std::string &filename, std::vector<UbTupleFloat3> &nodes,
                                         std::vector<UbTupleInt3> &triangles)
{
    string bobjFilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterBOBJ::writeTriangles to " << bobjFilename << " - start");

    gzFile gzf = gzopen(bobjFilename.c_str(), "wb1");

    size_t nofNodes     = nodes.size();
    size_t nofTriangles = triangles.size();

    // write to file
    size_t numVerts;
    // double v[3];
    if (sizeof(numVerts) != 4) {
        throw UbException(UB_EXARGS, "danger...");
    }
    numVerts = nofNodes;
    gzwrite(gzf, &numVerts, sizeof(numVerts));

    for (size_t k = 0; k < nofNodes; k++) {
        float vertp = val<1>(nodes[k]);
        gzwrite(gzf, &vertp, sizeof(vertp));
        vertp = val<2>(nodes[k]);
        gzwrite(gzf, &vertp, sizeof(vertp));
        vertp = val<3>(nodes[k]);
        gzwrite(gzf, &vertp, sizeof(vertp));
    }

    // NORMAL VECTOR
    // double n[3];
    gzwrite(gzf, &numVerts, sizeof(numVerts));
    for (size_t k = 0; k < nofNodes; k++) {
        // poly->GetPointData()->GetNormals()->GetTuple(k, n);
        float normp = 0.0; // n[0];
        gzwrite(gzf, &normp, sizeof(normp));
        normp = 0.0; // n[1];
        gzwrite(gzf, &normp, sizeof(normp));
        normp = 0.0; // n[2];
        gzwrite(gzf, &normp, sizeof(normp));
    }

    // vtkIdType npts = 3;
    // vtkIdType* pts;
    size_t numTris = nofTriangles;
    gzwrite(gzf, &numTris, sizeof(numTris));
    for (size_t k = 0; k < nofTriangles /*(size_t)poly->GetNumberOfPolys()*/; k++) {
        // poly->GetPolys()->GetNextCell(npts, pts);
        // int triIndex = *pts;
        // gzwrite(gzf, &triIndex, sizeof(triIndex));
        // triIndex = *(pts+1);
        // gzwrite(gzf, &triIndex, sizeof(triIndex));
        // triIndex = *(pts+2);
        // gzwrite(gzf, &triIndex, sizeof(triIndex));
        // poly->GetPolys()->GetNextCell(npts, pts);
        int triIndex = val<1>(triangles[k]); //*pts;
        gzwrite(gzf, &triIndex, sizeof(triIndex));
        triIndex = val<2>(triangles[k]); //*(pts+1);
        gzwrite(gzf, &triIndex, sizeof(triIndex));
        triIndex = val<3>(triangles[k]); //*(pts+2);
        gzwrite(gzf, &triIndex, sizeof(triIndex));
    }

    gzclose(gzf);

    UBLOG(logDEBUG1, "WbWriterBOBJ::writeTriangles to " << bobjFilename << " - end");

    return bobjFilename;
}
/*===============================================================================*/

#endif // CAB_ZLIB

//! \}
