#include "BoundaryValues.h"
#include <stdio.h>
#include <stdlib.h>


BoundaryValues::BoundaryValues(string ad)
{
	file.open(ad.c_str(), ios::in);

	if (!file) {
		cerr << "Fehler beim Oeffnen_Values" <<endl;
			exit(1);
	}
	init();
}


BoundaryValues::BoundaryValues(string ad, Parameter* para, string str)
{
	if (!file) {
		cerr << "Fehler beim Oeffnen_Values" <<endl;
		exit(1);
	}


	file.open(ad.c_str(), ios::in);
	if (!file) {
		para->setObj(str, false);
	} else {
		init();
		para->setObj(str, true);
	}
}


BoundaryValues::BoundaryValues(int neighbor, Parameter* para, string sor)
{
	string ad = para->getPossNeighborFiles(sor)[neighbor];
	file.open(ad.c_str(), ios::in);

	if (!file) {
		para->setIsNeighbor(false);
	} else {
		para->setIsNeighbor(true);
		init();
	}
}


BoundaryValues::BoundaryValues(int neighbor, Parameter* para, string sor, string dir)
{
	if (dir=="X")
	{
		string ad = para->getPossNeighborFilesX(sor)[neighbor];
		file.open(ad.c_str(), ios::in);

		if (file.fail()) {
			para->setIsNeighborX(false);
		} else {
			para->setIsNeighborX(true);
			init();
		}
	} 
	else if (dir=="Y")
	{
		string ad = para->getPossNeighborFilesY(sor)[neighbor];
		file.open(ad.c_str(), ios::in);

		if (file.fail()) {
			para->setIsNeighborY(false);
		} else {
			para->setIsNeighborY(true);
			init();
		}
	}
	else
	{
		string ad = para->getPossNeighborFilesZ(sor)[neighbor];
		file.open(ad.c_str(), ios::in);

		if (file.fail()) {
			para->setIsNeighborZ(false);
		} else {
			para->setIsNeighborZ(true);
			init();
		}
	}
}


BoundaryValues::~BoundaryValues(void)
{
	file.close();
}

void BoundaryValues::init() {

	string bufferString;
	unsigned int bufferInt;
	doubflo bufferDouble;


	getline(file,bufferString);
	bo=bufferString; // "bo" speichert die Art der Randbedingung

	getline(file,bufferString);//level einlesen
	level = atoi(bufferString.c_str());
	getline(file,bufferString);//ueberspringen der groesse

	int counter = 0;
	int laeuftmit = 0;
	while(counter <= 0 && laeuftmit <= level) { //schleife laeuft bis spaltenanzahl gefunden wurde
		getline(file, bufferString);
		for(unsigned int i = 0; i < bufferString.length(); i++)
		{
			if (bufferString[i] == ' ')
			{
				counter++;
			}
		}
		laeuftmit++;
	}
	int column = counter;
	//cout << "column: " << column << endl;


	if(bo=="noSlip"||bo=="slip") {

		vec1D_index.push_back(0);
		vec2D_index.push_back(vec1D_index);
		vec1D_index.clear();

		vec1D_column.push_back(0);
		vec2D_lvl.push_back(vec1D_column); 
		vec1D_column.clear(); 
		vec3D_data.push_back(vec2D_lvl);
		vec2D_lvl.clear();
	}

	else //if (bo == "velocity" || bo == "pressure" || bo == "periodic_x" || bo == "periodic_y" || bo == "periodic_z") {
	{
		file.clear();
		file.seekg (0, ios::beg);
		getline(file,bufferString); //Art wird uebersprungen
		getline(file,bufferString); //level wird uebersprungen

		////Levelgroessen aus der dritten Zeile werden eingelesen------------------------------------------------------------------//
		for(unsigned int i=0; i<=level; i++) {
			getline(file,bufferString);
			bufferInt = atoi(bufferString.c_str()); //eingelesene Zeile wird zum Integer gecastet
			vec_Size.push_back(bufferInt);
			if(bufferInt!=0) {//falls die Groesse Null ist, muss keine Zeile uebersprungen werden
				for(unsigned int j=0; j<bufferInt;j++){ 
					getline(file,bufferString);// ueberspringt alles bis zur naechsten Groesse
				}
			}
		}

		//Index wird in den Values eingelesen------------------------------------------------------------------------------------//
		file.clear();
		file.seekg (0, ios::beg);
		getline(file,bufferString); //Art wird uebersprungen
		getline(file,bufferString); //level wird uebersprungen

		vec3D_data.resize(level+1);
		for (unsigned int i = 0; i < vec3D_data.size(); i++) {
			vec3D_data[i].resize(column);
			for (unsigned int j = 0; j < vec3D_data[i].size(); j++) {
				vec3D_data[i][j].resize(vec_Size[i]);
			}
		}

		vec2D_index.resize(level+1);
		for (unsigned int i = 0; i < vec2D_index.size(); i++) {
			vec2D_index[i].resize(vec_Size[i]);
		}	

		//schleife zum Einlesen der Werte
		for (unsigned int j = 0; j <= level; j++) {
			if(vec_Size[j]==0) {
				getline(file,bufferString);
				continue;
			} //ist eine Groesse=0 -> continue

			getline(file,bufferString); // Groesse ueberspringen
			for (unsigned int elem = 0; elem < vec_Size[j]; elem++)
			{
				file >> bufferInt;
				vec2D_index[j][elem]=bufferInt;
				for (int col = 0; col < column; col++)
				{
					file >> bufferDouble;
					vec3D_data[j][col][elem] = bufferDouble;
				}
			}
			getline(file,bufferString);
		} // ende for-schleife zum Einlesen der Werte
	} // ende if
}


vector<unsigned int> BoundaryValues::getIndex(unsigned int level) {
	int zaehler = 0;
	for (vector<vector<unsigned int> >::iterator it = vec2D_index.begin() ; it != vec2D_index.end(); it++) {
		if ( zaehler==level) {return *it;}
		zaehler++;
	}
	cout<<"Levelgroesse nicht zulaessig"<<endl;
	exit(1);
}


vector<doubflo> BoundaryValues::getVec(unsigned int level, unsigned int column) {
	int zaehler1=0;
	int zaehler2=0;

	for (vector<vector<vector<doubflo> > >::iterator it = vec3D_data.begin() ; it != vec3D_data.end(); it++) {
		if (zaehler1==level) { 
			for(vector<vector<doubflo> >::iterator it2 = it->begin() ; it2 != it->end(); it2++) {
				if(zaehler2==column) {return *it2;}
				zaehler2++;
			}
		}	
		zaehler1++;
	}
	cout<<"Levelgroesse nicht zulaessig"<<endl;
	exit(1);
}


string BoundaryValues::getWay() {
	return this->bo;
}


void BoundaryValues::setProcNeighbor(bool pN){
	procNeighbor = pN;
}


bool BoundaryValues::getProcNeighbor(){
	return procNeighbor;
}
