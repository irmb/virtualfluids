#include "BoundaryQs.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;


BoundaryQs::BoundaryQs(string ad, bool binaer)
{

	if(binaer == true) {
		file.open(ad.c_str(), ios::in | ios::binary);

		if (!file) {
		cerr << "Fehler beim Oeffnen_Qs" <<endl;
			exit(1);
		}

		init_Binary();
	} else {
		file.open(ad.c_str(), ios::in );

		if (!file) {
		cerr << "Fehler beim Oeffnen_Qs" <<endl;
			exit(1);
		}

		init();
		
	}
}
BoundaryQs::BoundaryQs(string ad, Parameter* para, string str, bool binaer){

	if (binaer==true) file.open(ad.c_str(), ios::in | ios::binary);
	else file.open(ad.c_str(), ios::in);

	if (!file) {
		para->setObj(str, false);
	} else {
		if(binaer == true) {
			init_Binary();
		} else {
			init();
		}
		para->setObj(str, true);
	}
	
}
BoundaryQs::BoundaryQs(){
}



BoundaryQs::~BoundaryQs(void)
{
	file.close();
}

void BoundaryQs::init() {
	vector<uint32_t> vec1D_code;

	file.clear();
	file.seekg(0, ios::beg);
	string bufferString;
	unsigned int bufferInt;
	doubflo bufferDouble;


	getline(file, bufferString); //level einlesen
	level = atoi(bufferString.c_str());
	column = 27;
	vec_Size.resize(level + 1);
	vec3D_data.resize(level + 1);
	vec2D_index.resize(level + 1);


	//schleife zum Einlesen der Werte
	for (unsigned int j = 0; j <= level; j++) {
		vec3D_data[j].resize(column);
		getline(file, bufferString);
		bufferInt = atoi(bufferString.c_str()); //eingelesene Zeile wird zum Integer gecastet
		vec_Size[j] = bufferInt;
		for (int i = 0; i < column; i++) {
			vec3D_data[j][i].resize(vec_Size[j], -1);
		}
		vec2D_index[j].resize(vec_Size[j]);
		if (bufferInt == 0) {//falls die Groesse Null ist, muss keine Zeile uebersprungen werden
			continue;
		}
		vec1D_code.resize(vec_Size[j]);

		for (unsigned int elem = 0; elem < vec_Size[j]; elem++)
		{

			
			
			int zaehler = 26;
			file >> bufferInt;
			vec2D_index[j][elem] = bufferInt;
			//cout << "Index: " << bufferInt << endl;

			bufferInt = atoi(bufferString.c_str());
			file >> bufferInt;
			vec1D_code[elem] = bufferInt;
			//cout << vec1D_code[elem] << endl;
			//getline(file, bufferString); // Code ueberspringen

			while (vec1D_code[elem] != 0){

				if (vec1D_code[elem] % 2 == 1){

					file >> bufferDouble;

					vec3D_data[j][zaehler][elem] = bufferDouble;
				}
				vec1D_code[elem] /= 2;
				zaehler--;
			}
			getline(file, bufferString);

		}

		vec1D_code.clear();
	}// ende for-schleife zum Einlesen der Werte

}



void BoundaryQs::init_Binary() {

	vector<uint32_t> vec1D_code;

	file.clear();
	file.seekg(0, ios::beg);
	string bufferString;
	unsigned int bufferInt;
	doubflo bufferDouble;
	uint32_t bufferUint32_t;


	getline(file, bufferString); //level einlesen
	level = atoi(bufferString.c_str());
	column = 27;
	vec_Size.resize(level + 1);
	vec3D_data.resize(level + 1);
	vec2D_index.resize(level + 1);

	for (unsigned int j = 0; j <= level; j++) {
		vec3D_data[j].resize(column);
		getline(file, bufferString);
		bufferInt = atoi(bufferString.c_str()); 
		vec_Size[j] = bufferInt;
		for (int i = 0; i < column; i++) {
			vec3D_data[j][i].resize(vec_Size[j], -1);
		}
		vec2D_index[j].resize(vec_Size[j]);
		if (bufferInt == 0) {//falls die Groesse Null ist, muss keine Zeile uebersprungen werden
			continue;
		}
		vec1D_code.resize(vec_Size[j]);

		for (unsigned int elem = 0; elem < vec_Size[j]; elem++)
		{
		int zaehler = 26;
		file.read((char*)&bufferInt,sizeof(int));
		vec2D_index[j][elem] = bufferInt;
		

		file.read((char*)&bufferUint32_t,sizeof(uint32_t));
		vec1D_code[elem] = bufferUint32_t;
		while (vec1D_code[elem] != 0){

				if (vec1D_code[elem] % 2 == 1){

					file.read((char*)&bufferDouble,sizeof(double));
					vec3D_data[j][zaehler][elem] = bufferDouble;
				}
				vec1D_code[elem] /= 2;
				zaehler--;
		}
		getline(file, bufferString);
		}
		vec1D_code.clear();
	}
}

vector<vector<vector<doubflo> > > BoundaryQs::setgetBoundarys(vector<vector<vector<doubflo> > > qs) {

	int j=0;
	int i=0;
	for (vector<vector<vector<doubflo> > >::iterator it = vec3D_data.begin() ; it != vec3D_data.end(); it++) {
		i = 0;
		for(vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
			
			for(vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
				qs[j][i].push_back(*it3);
			}
			i++;
		}
		j++;
	}
	vec3D_data.clear();

	return qs;
}

vector<vector<vector<unsigned int> > > BoundaryQs::setgetBoundarys(vector<vector<vector<unsigned int> > > qs) {

	cout << "Laenge 1 vec3D_data: " << vec3D_data.size() << endl;
	cout << "Laenge 2 vec3D_data: " << vec3D_data[0].size() << endl;
	cout << "Laenge 3 vec3D_data: " << vec3D_data[0][0].size() << endl;	
	int j = 0;
	int i = 0;
	for (vector<vector<vector<doubflo> > >::iterator it = vec3D_data.begin(); it != vec3D_data.end(); it++) {
		i = 0;
		for (vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){

			for (vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
				qs[j][i].push_back(*it3);
			}
			i++;
		}
		j++;
	}
	vec3D_data.clear();

	return qs;
}




vector<vector<unsigned int> > BoundaryQs::setgetIndex(vector<vector<unsigned int> > index) {

	int i=0;

		for(vector<vector<unsigned int> >::iterator it2 = vec2D_index.begin(); it2 != vec2D_index.end(); it2++){

			for(vector<unsigned int>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
				if(*it3!=0) {
					index[i].push_back(*it3);
				}
			}
			i++;
		}

	vec2D_index.clear();
	return index;
}

void BoundaryQs::initArray(doubflo* ptr, unsigned int level, unsigned int column) {
	int zaehler=0;
	int zaehler2=0;
	int n=0;
	if (ptr != NULL) { //nur wenn genug Speicher vorhanden ist
		for (vector<vector<vector<doubflo> > >::iterator it = vec3D_data.begin() ; it != vec3D_data.end(); it++) {
			if(zaehler==level) {
				for(vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
					if(zaehler2==column) {
						for(vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
						ptr[n]=*it3;
						n++;
						}
					}
					zaehler2++; // zaehlt die Spalte mit
				}				
			}
			zaehler++; // zaehlt das Level mit
		}
		
	}
}

void BoundaryQs::initIndex(/*unsigned*/ int *ptr, unsigned int level) {
	int zaehler = 0;
	unsigned int n=0;
	if (ptr !=NULL) {
		for (vector<vector<unsigned int> >::iterator it = vec2D_index.begin() ; it != vec2D_index.end(); it++) {
			if ( zaehler==level) {
				for(vector<unsigned int>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
					ptr[n]=*it2;
					n++;
				}
			}
			zaehler++;
		}
	}
}

void BoundaryQs::initProcNeighbor(int* ptr, unsigned int level, unsigned int column) {
	int zaehler=0;
	int zaehler2=0;
	int n=0;
	if (ptr != NULL) { //nur wenn genug Speicher vorhanden ist
		for (vector<vector<vector<doubflo> > >::iterator it = vec3D_data.begin() ; it != vec3D_data.end(); it++) {
			if(zaehler==level) {
				for(vector<vector<doubflo> >::iterator it2 = it->begin(); it2 != it->end(); it2++){
					if(zaehler2==column) {
						for(vector<doubflo>::iterator it3 = it2->begin(); it3 != it2->end(); it3++){
							ptr[n]=(int)*it3;
							n++;
						}
					}
					zaehler2++; // zaehlt die Spalte mit
				}				
			}
			zaehler++; // zaehlt das Level mit
		}

	}
}

int BoundaryQs::getcolumn(){
	return this->column;
}

