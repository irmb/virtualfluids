#pragma once

#include <vector>
#include <fstream>

using namespace std;

class CoordNeighborGeoV 
{
protected:
	ifstream file; 
	unsigned int level; 
	vector<unsigned int> vec_Size; // Die Groeßen der Felder

	vector<unsigned int> vec_data; //das Feld mit Werten, temporär zum Füllen von vec2D_data
	vector<doubflo> vec_data_Coord; //das Feld mit Werten, temporär zum Füllen von vec2D_data
	vector< vector<unsigned int> > vec2D_data; //alle Felder mit Werten gegliedert nach Level
	vector< vector< doubflo> > vec2D_data_Coord; //alle Felder mit Werten gegliedert nach Level

public:
	CoordNeighborGeoV();
	CoordNeighborGeoV(string ad, bool binaer, bool coord);
	~CoordNeighborGeoV(void);
	void init(bool coord); //füllt die Vektoren 
	void init_Binary(bool coord);
	

	unsigned int getLevel();//liefert einen Wert, die größte Levelzahl
	unsigned int getSize(unsigned int level); //liefert die Größe des gesuchten Levels
	vector<unsigned int>getVec(unsigned int level);
	void setVec(unsigned int level, vector<unsigned int> vec);

	void initArray(unsigned int *int_ptr, unsigned int level );
	void initArrayCoord(doubflo *int_ptr, unsigned int level );
};

