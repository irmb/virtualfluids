#include "OffsetScale.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

OffsetScale::OffsetScale(string ad, bool off)
{
	file.open(ad.c_str(), ios::in);

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
	string buffer;
	unsigned int bufferInt;

    this->readLevel();

	//Schleife zum Einlesen der Levelgroessen
	for(unsigned int i=1; i<= maxLevel; i++) {
		getline(file,buffer);
		unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
        levelSizes.push_back(bufferInt);
		getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
	}
	
	file.clear();
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	getline(file,buffer); //level wird ignoriert


	//einlesen der Werte
    scale.resize(maxLevel +1);
	for(unsigned lvl = 0; lvl < maxLevel; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		getline(file,buffer); // Groesse ignorieren
		for (unsigned int i = 0; i < levelSizes[lvl]; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		{
			file >> bufferInt;
            scale[lvl].push_back(bufferInt);
		}
		getline(file, buffer);
	}
}

void OffsetScale::initOffset() {
	file.clear();
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	//Level aus der ersten Zeile wird ausgelesen
	string buffer;
	real bufferDouble;

    this->readLevel();

	//Schleife zum Einlesen der Levelgroessen
	for(unsigned int i=1; i<= maxLevel; i++) {
		getline(file,buffer);
		unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
        levelSizes.push_back(bufferInt);
		getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
	}
	
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	getline(file,buffer); //level wird ignoriert

	//einlesen der werte
	vec2D_dataOffset.resize(maxLevel +1);
	for(unsigned lvl = 0; lvl < maxLevel; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		getline(file,buffer); // Groesse ignorieren
		for (unsigned int i = 0; i < levelSizes[lvl]*3; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
														  /////////////Unterschied zu Scale:  vec_Size[lvl]*3       //////////////////////
		{
			file >> bufferDouble;
			vec2D_dataOffset[lvl].push_back(bufferDouble);
		}
		getline(file, buffer);
	}
}


void OffsetScale::initArrayOffset(real *x_ptr,real *y_ptr,real *z_ptr, unsigned int level) {
	
	int zaehler=0;
	int n=0;
	int x_help=0;
	int y_help=1;
	int z_help=2;
	int x=0;
	int y=0;
	int z=0;

	if (x_ptr != NULL && y_ptr !=NULL && z_ptr !=NULL) {
		for (vector<vector<real> >::iterator it = vec2D_dataOffset.begin() ; it != vec2D_dataOffset.end(); it++) {		
			if(zaehler==level) {
				for(vector<real>::iterator it2 = it->begin(); it2 != it->end(); it2++){
					if(n==x_help){
						x_ptr[x]=*it2;
						x_help=x_help+3;
						x++;
					} else if(n==y_help){
						y_ptr[y]=*it2;
						y_help=y_help+3;
						y++;
					} else if(n==z_help){
						z_ptr[z]=*it2;
						z_help=z_help+3;
						z++;
					}
					n++;
				}
			}
			zaehler++;	
		}

	}

}

void OffsetScale::initScale(unsigned int* data, unsigned int level)
{
    for (int index = 0; index < scale[level].size(); index++)
        data[index] = scale[level][index];
}

