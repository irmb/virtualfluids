#include "OffsetScale.h"
#include <stdio.h>
#include <stdlib.h>

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

	getline(file,buffer);
	level = atoi(buffer.c_str());

	//Schleife zum Einlesen der Levelgroessen
	for(unsigned int i=1; i<=level; i++) {
		getline(file,buffer);
		unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
		vec_Size.push_back(bufferInt);
		getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
	}
	
	file.clear();
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	getline(file,buffer); //level wird ignoriert


	//einlesen der Werte
	vec2D_data.resize(level+1);
	for(unsigned lvl = 0; lvl < level; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		getline(file,buffer); // Groesse ignorieren
		for (unsigned int i = 0; i < vec_Size[lvl]; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		{
			file >> bufferInt;
			vec2D_data[lvl].push_back(bufferInt);
		}
		getline(file, buffer);
	}
}

void OffsetScale::initOffset() {
	file.clear();
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	//Level aus der ersten Zeile wird ausgelesen
	string buffer;
	doubflo bufferDouble;

	getline(file,buffer);
	stringstream s1(buffer);
	s1 >> level;

	//Schleife zum Einlesen der Levelgroessen
	for(unsigned int i=1; i<=level; i++) {
		getline(file,buffer);
		unsigned int bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
		vec_Size.push_back(bufferInt);
		getline(file,buffer); //die Zeile mit den Koordinaten muss uebersprungen werden
	}
	
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	getline(file,buffer); //level wird ignoriert

	//einlesen der werte
	vec2D_dataOffset.resize(level+1);
	for(unsigned lvl = 0; lvl < level; lvl++){/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
		getline(file,buffer); // Groesse ignorieren
		for (unsigned int i = 0; i < vec_Size[lvl]*3; i++)/////////////Unterschied zu CoordNeighborGeoV:  < statt <=//////////////////////
														  /////////////Unterschied zu Scale:  vec_Size[lvl]*3       //////////////////////
		{
			file >> bufferDouble;
			vec2D_dataOffset[lvl].push_back(bufferDouble);
		}
		getline(file, buffer);
	}
}


void OffsetScale::initArrayOffset(doubflo *x_ptr,doubflo *y_ptr,doubflo *z_ptr, unsigned int level) {
	
	int zaehler=0;
	int n=0;
	int x_help=0;
	int y_help=1;
	int z_help=2;
	int x=0;
	int y=0;
	int z=0;

	if (x_ptr != NULL && y_ptr !=NULL && z_ptr !=NULL) {
		for (vector<vector<doubflo> >::iterator it = vec2D_dataOffset.begin() ; it != vec2D_dataOffset.end(); it++) {		
			if(zaehler==level) {
				for(vector<doubflo>::iterator it2 = it->begin(); it2 != it->end(); it2++){
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