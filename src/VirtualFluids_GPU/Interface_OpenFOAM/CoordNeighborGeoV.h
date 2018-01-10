#pragma once

#include <vector>
#include <fstream>

using namespace std;

class CoordNeighborGeoV 
{
protected:
	ifstream file; 
	unsigned int level; 
	vector<unsigned int> vec_Size; // Die Groe�en der Felder

	vector<unsigned int> vec_data; //das Feld mit Werten, tempor�r zum F�llen von vec2D_data
	vector<doubflo> vec_data_Coord; //das Feld mit Werten, tempor�r zum F�llen von vec2D_data
	vector< vector<unsigned int> > vec2D_data; //alle Felder mit Werten gegliedert nach Level
	vector< vector< doubflo> > vec2D_data_Coord; //alle Felder mit Werten gegliedert nach Level

public:
	CoordNeighborGeoV();
	CoordNeighborGeoV(string ad, bool binaer, bool coord);
	~CoordNeighborGeoV(void);
	void init(bool coord); //f�llt die Vektoren 
	void init_Binary(bool coord);
	

	unsigned int getLevel();//liefert einen Wert, die gr��te Levelzahl
	unsigned int getSize(unsigned int level); //liefert die Gr��e des gesuchten Levels
	vector<unsigned int>getVec(unsigned int level);
	void setVec(unsigned int level, vector<unsigned int> vec);

	void initArray(unsigned int *int_ptr, unsigned int level );
	void initArrayCoord(doubflo *int_ptr, unsigned int level );
};

