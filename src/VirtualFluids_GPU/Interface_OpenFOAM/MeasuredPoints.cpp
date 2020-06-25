#include "MeasuredPoints.h"
#include <stdlib.h>

MeasuredPoints::MeasuredPoints(void)
{
}

MeasuredPoints::MeasuredPoints(string ad){
	file.open(ad.c_str(), ios::in | ios::binary);

	if (!file) {
		cerr << "Fehler beim Oeffnen Measured Points" <<endl;
			exit(1);
	}

	this->init();		

}

MeasuredPoints::~MeasuredPoints(void)
{
}




void MeasuredPoints::init() {
	
	string bufferString;
	unsigned int bufferInt;

	getline(file,bufferString);
	getline(file,bufferString);//level einlesen
	level = atoi(bufferString.c_str());

	this->vec_Size.resize(level);
	this->vec2D_data.resize(level);

	for (int i=0; i<level;i++) {
		getline(file,bufferString);
		bufferInt = atoi(bufferString.c_str()); 

		this->vec_Size[i]=bufferInt;

		this->vec2D_data[i].resize(vec_Size[i]);
		if(vec_Size[i] != 0) {
			for ( int j=0; j<vec_Size[i]; j++) {
				getline(file,bufferString);
				bufferInt = atoi(bufferString.c_str()); 
				this->vec2D_data[i][j]=bufferInt;
			}
		}


	}



}