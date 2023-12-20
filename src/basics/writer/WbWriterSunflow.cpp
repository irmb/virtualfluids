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
#include <basics/utilities/UbLogger.h>
#include <basics/writer/WbWriterSunflow.h>
#include <cstring>

using namespace std;

/*===============================================================================*/
std::string WbWriterSunflow::writeTriangles(const string &filename, vector<UbTupleFloat3> &nodes,
                                            vector<UbTupleInt3> &triangles)
{
    string sunflowFilename = filename + getFileExtension();
    UBLOG(logDEBUG1, "WbWriterSunflow::writeTriangles to " << sunflowFilename << " - start");

    std::ofstream out(sunflowFilename.c_str());
    if (!out) {
        out.clear(); // flags ruecksetzen (ansonsten liefert utern if(!out) weiterhin true!!!
        string path = UbSystem::getPathFromString(sunflowFilename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            out.open(sunflowFilename.c_str());
        }
        if (!out)
            throw UbException(UB_EXARGS, "couldn't open file " + sunflowFilename);
    }

    // General part

    // Image details
    out << "image {" << endl;
    out << "   resolution 640 480" << endl;
    out << "   aa 0 1" << endl;
    out << "   filter mitchell" << endl;
    out << "}" << endl << endl;

    // Camera position
    out << "camera {" << endl;
    out << "   type pinhole" << endl;
    out << "   eye    -0.25 -0.3 0.13" << endl;
    out << "   target -0.1 0.1 0.13" << endl;
    out << "   up     0 0 1" << endl;
    out << "   fov    60" << endl;
    out << "   aspect 1.333333" << endl;
    out << "}" << endl << endl;

    // Light
    out << "light {" << endl;
    out << "   type ibl" << endl;
    out << "   image sky_small.hdr" << endl;
    out << "   center 0 -1 0" << endl;
    out << "   up 0 0 1" << endl;
    out << "   lock true" << endl;
    out << "   samples 200" << endl;
    out << "}" << endl << endl;

    // Shaders
    out << "shader {" << endl;
    out << "   name default-shader" << endl;
    out << "   type diffuse" << endl;
    out << "   diff 0.25 0.25 0.25" << endl;
    out << "}" << endl << endl;

    out << "shader {" << endl;
    out << "   name Glass" << endl;
    out << "   type glass" << endl;
    out << "   eta 1.333" << endl;
    out << "   color 0.1 0.3 0.8" << endl;
    out << "}" << endl << endl;

    out << "shader {" << endl;
    out << "   name Mirror" << endl;
    out << "   type mirror" << endl;
    out << "   refl 0.7 0.7 0.7" << endl;
    out << "}" << endl << endl;

    // Objects
    // a) Ground plane
    out << "object {" << endl;
    out << "   shader default-shader" << endl;
    out << "   type plane" << endl;
    out << "   p 0 0 0" << endl;
    out << "   n 0 0 1" << endl;
    out << "}" << endl << endl;

    // b) Mesh
    out << "object {" << endl;
    out << "   shader Glass" << endl;
    out << "   transform {" << endl;
    out << "      rotatey 270.0" << endl;
    out << "   }" << endl;
    out << "   type generic-mesh" << endl;
    out << "      name polySurfac" << endl << endl;

    // POINTS SECTION
    int nofNodes = (int)nodes.size();
    out << "   points " << nofNodes << endl;
    for (int n = 0; n < nofNodes; n++)
        out << "      " << val<1>(nodes[n]) << " " << val<2>(nodes[n]) << " " << val<3>(nodes[n]) << endl;

    // TRIANGLES SECTION
    int nofTriangles = (int)triangles.size();
    out << "   triangles " << nofTriangles << endl;
    for (int c = 0; c < nofTriangles; c++)
        out << "      " << val<1>(triangles[c]) << " " << val<2>(triangles[c]) << " " << val<3>(triangles[c]) << endl;

    // FOOTER
    out << "   normals none" << endl;
    out << "   uvs none" << endl;
    out << "}" << endl;

    out.close();
    UBLOG(logDEBUG1, "WbWriterSunflow::writeTriangles to " << sunflowFilename << " - end");

    return sunflowFilename;
}
/*===============================================================================*/

//! \}
