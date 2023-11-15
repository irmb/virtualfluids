#include "OffsetScale.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

OffsetScale::OffsetScale(std::string ad, bool off)
{
    file.open(ad.c_str(), std::ios::in);

    if (!file) {
        cerr << "Fehler beim Oeffnen" <<endl;
            exit(1);
    }
    if(off==true){
        initOffset();
    }else {
        init();
    }
}
OffsetScale::~OffsetScale(void)
{
    file.close();
}

void OffsetScale::init() {
    //Level aus der ersten Zeile wird ausgelesen
    //string buffer;
    //unsigned int bufferInt;

    this->readLevel();
    levelSizes.resize(maxLevel);
    scale.resize(maxLevel);

    for (uint level = 0; level < maxLevel; level++)
    {
        readLevelSize(level);
        scale[level].resize(levelSizes[level]);
        for (uint index = 0; index < levelSizes[level]; index++)
            file >> scale[level][index];
    }

    ////Schleife zum Einlesen der Levelgroessen
    //for(unsigned int i=1; i<= maxLevel; i++) {
    //    getline(file,buffer);
    //    unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
 //       levelSizes.push_back(bufferInt);
    //    getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
    //}
    //
    //file.clear();
    //file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
    //getline(file,buffer); //level wird ignoriert


    ////einlesen der Werte
 //   scale.resize(maxLevel +1);
    //for(unsigned lvl = 0; lvl < maxLevel; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
    //    getline(file,buffer); // Groesse ignorieren
    //    for (unsigned int i = 0; i < levelSizes[lvl]; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
    //    {
    //        file >> bufferInt;
 //           scale[lvl].push_back(bufferInt);
    //    }
    //    getline(file, buffer);
    //}
}

void OffsetScale::initOffset()
{

    this->readLevel();
    levelSizes.resize(maxLevel);
    offset.resize(maxLevel);

    for (uint level = 0; level < maxLevel; level++)
    {
        readLevelSize(level);
        offset[level].resize(levelSizes[level] * 3);
        for (uint index = 0; index < levelSizes[level] * 3; index++)
            file >> offset[level][index];
    }





    //file.clear();
    //file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
    ////Level aus der ersten Zeile wird ausgelesen
    //string buffer;
    //real bufferDouble;

 //   this->readLevel();

    ////Schleife zum Einlesen der Levelgroessen
    //for(unsigned int i=1; i<= maxLevel; i++) {
    //    getline(file,buffer);
    //    unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
 //       levelSizes.push_back(bufferInt);
    //    getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
    //}
    //
    //file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
    //getline(file,buffer); //level wird ignoriert

    ////einlesen der werte
    //offset.resize(maxLevel +1);
    //for(unsigned lvl = 0; lvl < maxLevel; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
    //    getline(file,buffer); // Groesse ignorieren
    //    for (unsigned int i = 0; i < levelSizes[lvl]*3; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
    //                                                      /////////////Unterschied zu Scale:  vec_Size[lvl]*3       //////////////////////
    //    {
    //        file >> bufferDouble;
    //        offset[lvl].push_back(bufferDouble);
    //    }
    //    getline(file, buffer);
    //}
}


void OffsetScale::initArrayOffset(real *x_ptr,real *y_ptr,real *z_ptr, uint level)
{
    int coordIndex = 0;
    for (std::size_t index = 0; index < offset[level].size(); index+=3)
    {
        x_ptr[coordIndex] = offset[level][index];
        y_ptr[coordIndex] = offset[level][index + 1];
        z_ptr[coordIndex] = offset[level][index + 2];
        coordIndex++;
    }
}

void OffsetScale::initScale(unsigned int* data, unsigned int level)
{
    for (std::size_t index = 0; index < scale[level].size(); index++)
        data[index] = scale[level][index];
}

