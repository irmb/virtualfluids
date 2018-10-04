#define _CRT_SECURE_NO_DEPRECATE
#include "STLReader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#include <GridGenerator/geometries/Vertex/Vertex.h>
#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>

#include <utilities/logger/Logger.h>


std::vector<Triangle> STLReader::readSTL(const std::string& name)
{
    std::ifstream file(name.c_str());
    if (file.is_open()) {
        std::string line;
        std::getline(file, line);
        line[strcspn(line.c_str(), "\r\n")] = 0;
        if (strcmp(line.c_str(), "solid ascii") == 0) {
            file.close();
            *logging::out << logging::Logger::INFO_INTERMEDIATE << "start reading ascii STL file: " + name + "\n";
            return readASCIISTL(name);
        }
        else {
            file.close();
            *logging::out << logging::Logger::INFO_INTERMEDIATE << "start reading binary STL file: " + name + "\n";
            return readBinarySTL(name);
        }
    }

     *logging::out << logging::Logger::INFO_INTERMEDIATE << "can't open STL-file" + name + " ... exit program! \n";
     exit(1);
}

std::vector<Triangle> STLReader::readSTL(const std::string & name, FileType fileType, const std::vector<uint> ignorePatches)
{
    if ( fileType == ascii ) return readASCIISTLWithPatches(name, ignorePatches);
    else                     return readBinarySTL(name);
}


std::vector<Triangle> STLReader::readASCIISTL(const std::string& name)
{
    const int lines = countLinesInFile(name);
    const int nTriangles = (lines) / 7; // seven lines per triangle

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Number of Triangles: " << nTriangles << "\n";
    std::vector<Triangle> triangles;

    std::string line;
    std::ifstream file;
    file.open(name.c_str(), std::ifstream::in);
    std::getline(file, line); // solid ascii

    for (int t = 0; t < nTriangles; t++) {
        Vertex normal = parseLineToCoordinates(file, "%*s %*s %f %f %f");
        getline(file, line); // outer loop
        Vertex p1 = parseLineToCoordinates(file, "%*s %f %f %f");
        Vertex p2 = parseLineToCoordinates(file, "%*s %f %f %f");
        Vertex p3 = parseLineToCoordinates(file, "%*s %f %f %f");
        getline(file, line); //endloop
        getline(file, line); //endfacet

        Triangle tri = Triangle(p1, p2, p3, normal);
        tri.calcNormal();
        triangles.push_back(tri);
    }
    file.close();
    return triangles;
}


std::vector<Triangle> STLReader::readASCIISTLWithPatches(const std::string& name, const std::vector<uint> ignorePatches)
{
    *logging::out << logging::Logger::INFO_HIGH << "Start reading ascii STL file:\n";
    *logging::out << logging::Logger::INFO_HIGH << "    " + name + "\n";

    std::vector<Triangle> triangles;

    std::string line;
    std::ifstream file;
    file.open(name.c_str(), std::ifstream::in);

    uint currentPatchIndex = 0;

    uint currentFacetLine = 0;

    bool ignoreCurrentPatch = false;

    Vertex vertex1, vertex2, vertex3, normal;

    while( std::getline(file, line) ){

        // trim the string
        line = line.substr( line.find_first_not_of(" "), line.find_last_not_of(" ") + 1 );
        
        if( line.substr( 0, line.find(' ') ) == "color" ) continue;

        // ========================================================================================
        if     ( currentFacetLine == 0 && line.substr( 0, line.find(' ') ) == "solid" )
        {
            ignoreCurrentPatch = std::find( ignorePatches.begin(), ignorePatches.end(), currentPatchIndex ) != ignorePatches.end();

            if( !ignoreCurrentPatch )
                *logging::out << logging::Logger::INFO_INTERMEDIATE << "    Reading STL-Group " << line.substr( line.find(' ') + 1 ) << " as patch " << currentPatchIndex << "\n";
            else
                *logging::out << logging::Logger::WARNING           << "    Ignoring STL-Group " << line.substr( line.find(' ') + 1 ) << " as patch " << currentPatchIndex << "\n";

            currentFacetLine++;
        }
        else if( currentFacetLine == 1 && line.substr( 0, line.find(' ') ) == "endsolid" )
        {
            currentFacetLine = 0;
            currentPatchIndex++;
        }
        // ========================================================================================
        else if( currentFacetLine == 1 && line.substr( 0, line.find(' ') ) == "facet" )
        {
            normal = parseLineToCoordinates(line, "%*s %*s %f %f %f");
            currentFacetLine++;
        }
        else if( currentFacetLine == 2 && line.substr( 0, line.find(' ') ) == "outer" )
        {
            currentFacetLine++;
        }
        else if( currentFacetLine == 3 && line.substr( 0, line.find(' ') ) == "vertex" )
        {
            vertex1 = parseLineToCoordinates(line, "%*s %f %f %f");
            currentFacetLine++;
        }
        else if( currentFacetLine == 4 && line.substr( 0, line.find(' ') ) == "vertex" )
        {
            vertex2 = parseLineToCoordinates(line, "%*s %f %f %f");
            currentFacetLine++;
        }
        else if( currentFacetLine == 5 && line.substr( 0, line.find(' ') ) == "vertex" )
        {
            vertex3 = parseLineToCoordinates(line, "%*s %f %f %f");
            currentFacetLine++;
        }
        else if( currentFacetLine == 6 && line.substr( 0, line.find(' ') ) == "endloop" )
        {
            currentFacetLine++;
        }
        else if( currentFacetLine == 7 && line.substr( 0, line.find(' ') ) == "endfacet" )
        {
            if( !ignoreCurrentPatch ){
                Triangle tri = Triangle(vertex1, vertex2, vertex3, normal);
                tri.calcNormal();

                tri.patchIndex = currentPatchIndex;

                triangles.push_back(tri);
            }

            currentFacetLine = 1;
        }
        else
        {
            throw std::runtime_error("STL-File does not comply with standard!\n");
        }
    }

    file.close();

    *logging::out << logging::Logger::INFO_HIGH << "Done reading ascii STL file\n";

    return triangles;
}

std::vector<Triangle> STLReader::readBinarySTL(const std::string& name)
{
    const std::string mode = "rb";
    FILE *file = fopen(name.c_str(), mode.c_str());

    char header_info[80] = "";
    fread(header_info, sizeof(char), 80, file);

    char nTri[4];
    fread(nTri, sizeof(char), 4, file);
    unsigned long nTriLong = *((unsigned long*)nTri);
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Number of Triangles: " << nTriLong << "\n";
    std::vector<Triangle> triangles;

    char facet[50];
    for (unsigned int i = 0; i < nTriLong; i++){
        fread(facet, sizeof(char), 50, file);

        Vertex normal = getVertexFromChar(facet);

        Vertex p1 = getVertexFromChar(facet + 12);
        Vertex p2 = getVertexFromChar(facet + 24);
        Vertex p3 = getVertexFromChar(facet + 36);

        triangles.push_back(Triangle(p1, p2, p3, normal));
    }
	fclose(file);

    return triangles;
}

std::vector<Triangle> STLReader::readSTL(const BoundingBox &box, const std::string& name)
{
	std::ifstream file(name.c_str());
	if (file.is_open()) {
		std::string line;
		std::getline(file, line);
		line[strcspn(line.c_str(), "\r\n")] = 0;
		if (strcmp(line.c_str(), "solid ascii") == 0) {
			file.close();
			*logging::out << logging::Logger::INFO_INTERMEDIATE << "start reading ascii STL file: " + name + "\n";
			return readASCIISTL(box, name);
		}
		else {
			file.close();
			*logging::out << logging::Logger::INFO_INTERMEDIATE << "start reading binary STL file: " + name + "\n";
			std::vector<Triangle> triangles = readBinarySTL(box, name);
			return triangles;
		}
	}
	else {
		*logging::out << logging::Logger::INFO_INTERMEDIATE << "can't open STL-file" + name + "\n";
		exit(1);
	}
}

std::vector<Triangle> STLReader::readASCIISTL(const BoundingBox &box, const std::string& name)
{
    int lines = countLinesInFile(name);
    int nTriangles = (lines) / 7; // seven lines per triangle
    std::cout << "Number of Triangles: " << nTriangles << std::endl;
    std::vector<Triangle> triangles;

    std::string line;
    std::ifstream file;
    file.open(name.c_str(), std::ifstream::in);
    std::getline(file, line); // solid ascii

    for (int i= 0; i < nTriangles; i++) {
        Vertex normal = parseLineToCoordinates(file, "%*s %*s %f %f %f");
        getline(file, line); // outer loop
        Vertex p1 = parseLineToCoordinates(file, "%*s %f %f %f");
        Vertex p2 = parseLineToCoordinates(file, "%*s %f %f %f");
        Vertex p3 = parseLineToCoordinates(file, "%*s %f %f %f");
        getline(file, line); //endloop
        getline(file, line); //endfacet

        Triangle t(p1, p2, p3, normal);
        t.calcNormal();
        if (box.isInside(t) || box.intersect(t))
            triangles.push_back(t);
    }
    file.close();
    return triangles;
}


std::vector<Triangle> STLReader::readBinarySTL(const BoundingBox &box, const std::string& name)
{
    FILE *file;
    std::string mode = "rb";
    file = fopen(name.c_str(), mode.c_str());

    char header_info[80] = "";
    char nTri[4];
    unsigned long nTriLong;
  
    fread(header_info, sizeof(char), 80, file);


    fread(nTri, sizeof(char), 4, file);
    nTriLong = *((unsigned long*)nTri);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Number of Triangles complete geometry: " << nTriLong << "\n";
    std::vector<Triangle> triangles;

    char facet[50];
    for (unsigned int i = 0; i < nTriLong; i++){
        fread(facet, sizeof(char), 50, file);

        Vertex normal = getVertexFromChar(facet);

        Vertex p1 = getVertexFromChar(facet + 12);
        Vertex p2 = getVertexFromChar(facet + 24);
        Vertex p3 = getVertexFromChar(facet + 36);

        Triangle t(p1, p2, p3, normal);
        if (box.isInside(t) || box.intersect(t))
            triangles.push_back(t);
    }
    int size = triangles.size();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Number of Triangles in process: " << size << "\n";
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Complete reading STL file. \n";

	fclose(file);

    return triangles;
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
Vertex STLReader::parseLineToCoordinates(std::ifstream& file, std::string format)
{
    std::string line;
    getline(file, line);
    const char* buffer = line.c_str();
    float x, y, z;
    sscanf(buffer, format.c_str(), &x, &y, &z);
    return Vertex(x, y, z);
}

Vertex STLReader::parseLineToCoordinates(const std::string& line, const std::string format)
{
    const char* buffer = line.c_str();
    float x, y, z;
    sscanf(buffer, format.c_str(), &x, &y, &z);
    return Vertex(x, y, z);
}

int STLReader::countLinesInFile(std::string name)
{
    std::ifstream file;
    file.open(name.c_str(), std::ifstream::in);
    int nTriLong = (int)std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
    file.close();
    return nTriLong;
}

Vertex STLReader::getVertexFromChar(const char* facet)
{
    char f1[4] = { facet[0],
        facet[1], facet[2], facet[3] };

    char f2[4] = { facet[4],
        facet[5], facet[6], facet[7] };

    char f3[4] = { facet[8],
        facet[9], facet[10], facet[11] };

    float xx = *((float*)f1);
    float yy = *((float*)f2);
    float zz = *((float*)f3);

    return Vertex((real)(xx), (real)(yy), (real)(zz));
}
