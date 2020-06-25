#include "CoordNeighborGeoV.h"
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

CoordNeighborGeoV::CoordNeighborGeoV(void){
}

CoordNeighborGeoV::CoordNeighborGeoV(string ad, bool binaer, bool coord){
	file.open(ad.c_str(), ios::in | ios::binary);

	if (!file) {
		cerr << "Fehler beim Oeffnen_CoordNeighborGeo" <<endl;
			exit(1);
	}
	if(binaer==true) {
		init_Binary(coord);
	} else {
		init(coord);
	}
		
}

CoordNeighborGeoV::~CoordNeighborGeoV(){
	file.close();
}

void CoordNeighborGeoV::init(bool coord) {
	//Level aus der ersten Zeile wird eingelesen
	string buffer;
	unsigned int bufferInt;
	doubflo bufferDoubflo;

	getline(file,buffer);
	level = atoi(buffer.c_str());

	//Schleife zum Einlesen der Levelgroessen
	for(unsigned int i=0; i<=level; i++) {
		getline(file,buffer);
		bufferInt = atoi(buffer.c_str()); //eingelesene Zeile wird zum Integer gecastet
		vec_Size.push_back(bufferInt);
		getline(file,buffer);
	}
	file.clear();
	file.seekg (0, ios::beg); // file wird wieder auf den Anfang gesetzt
	getline(file,buffer); //level wird ignoriert

	
	//einlesen der Werte
	if(coord == true) {
		vec2D_data_Coord.resize(level+1);
		for(unsigned lvl = 0; lvl <= level; lvl++){
			getline(file,buffer); // Groesse ignorieren
			for (unsigned int i = 0; i <= vec_Size[lvl]; i++)
			{
				file >> bufferDoubflo;
				vec2D_data_Coord[lvl].push_back(bufferDoubflo);
			}
			getline(file, buffer);
		}
	} else {
		vec2D_data.resize(level+1);
		for(unsigned lvl = 0; lvl <= level; lvl++){
			getline(file,buffer); // Groesse ignorieren
			for (unsigned int i = 0; i <= vec_Size[lvl]; i++)
			{
				file >> bufferInt;
				vec2D_data[lvl].push_back(bufferInt);
			}
			getline(file, buffer);
		}
	} // ende else
} // ende methode


void CoordNeighborGeoV::init_Binary(bool coord) {
	char c;
	string buffer;
	unsigned int bufferInt;
	double bufferDoubflo;

	//level einlesen
	getline(file,buffer);
	level = atoi(buffer.c_str());
	//Groesse des ersten Levels einlesen
	getline(file,buffer);
	bufferInt = atoi(buffer.c_str());
	vec_Size.push_back(bufferInt);
	
	
	for(unsigned int i=0; i<=level;i++) {
		if(coord == true) {
			file >> bufferDoubflo;
			vec_data_Coord.push_back(bufferDoubflo);//erste null
			file.get(c);//ueberspringt leerzeichen
		
			for(unsigned int j=0; j<vec_Size[i]; j++) {
				file.read((char*)&bufferDoubflo,sizeof(double));
				vec_data_Coord.push_back((doubflo)bufferDoubflo);
			}

			vec2D_data_Coord.push_back(vec_data_Coord);
			vec_data_Coord.clear();
		} else {
			file >> bufferInt;
			vec_data.push_back(bufferInt);//erste null
			file.get(c);//ueberspringt leerzeichen
		
			for(unsigned int j=0; j<vec_Size[i]; j++) {
				file.read((char*)&bufferInt,sizeof(unsigned int));
				vec_data.push_back(bufferInt);
			}

			vec2D_data.push_back(vec_data);
			vec_data.clear();
		}

		if(i==level) break;
		file.get(c);//ueberspringt Leerzeichen
		getline(file,buffer);
		bufferInt = atoi(buffer.c_str());
		vec_Size.push_back(bufferInt);
	}
}


unsigned int CoordNeighborGeoV::getLevel() {
	return level;
}

unsigned int CoordNeighborGeoV::getSize(unsigned int level) {
	int zaehler=0;
	for (vector<unsigned int>::iterator it = vec_Size.begin() ; it != vec_Size.end(); it++)	
	{if (zaehler == level){return *it; break;}// durchlauft den Vektor mit den Größen. Break wenn level dem Schleifendurchgang entspricht
	zaehler++;
	
	}
	 //cout<<"Levelgroesse nicht zulaessig"<<endl;
	 cout<<"Levelgroesse nicht zulaessig zaehler: "<< zaehler << " level: "<< level << " vec_Size.size(): " << vec_Size.size() <<endl;
	 exit(1);
}


void CoordNeighborGeoV::initArray(unsigned int *int_ptr, unsigned int level) {

	int zaehler=0;
	unsigned int n=0;
	if (int_ptr != NULL) { //nur wenn genug Speicher vorhanden ist
		for (vector<vector<unsigned int> >::iterator it = vec2D_data.begin() ; it != vec2D_data.end(); it++) {
			if(zaehler==level) {
				for(vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++){
					int_ptr[n]=*it2;
						n++;
				}
			}
			zaehler++;	
		}
	}
}

vector<unsigned int> CoordNeighborGeoV::getVec(unsigned int level) {
	int zaehler=0;
	for (vector<vector<unsigned int> >::iterator it = vec2D_data.begin() ; it != vec2D_data.end(); it++) {
		if (zaehler==level) {return *it;}
		zaehler++;	
	}
	cout<<"Levelgroesse nicht zulaessig"<<endl;
	 exit(1);
}

void CoordNeighborGeoV::setVec(unsigned int level, vector<unsigned int> vec) {
	this->vec2D_data[level]=vec;
	//for (int i=0; i<=2200; i++) {
		//cout <<"Test im Setter: "<< i <<": " << vec[i] << endl;
	//}
}







void CoordNeighborGeoV::initArrayCoord(doubflo *int_ptr, unsigned int level) {

	int zaehler=0;
	unsigned int n=0;
	if (int_ptr != NULL) { //nur wenn genug Speicher vorhanden ist
		for (vector<vector<doubflo> >::iterator it = vec2D_data_Coord.begin() ; it != vec2D_data_Coord.end(); it++) {		
			if(zaehler==level) {
				for(vector<doubflo>::iterator it2 = it->begin(); it2 != it->end(); it2++){
					int_ptr[n]=*it2;
					n++;
				}
			}
			zaehler++;	
		}
	}
}
