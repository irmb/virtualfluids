#include "MeasuredPoints.h"
#include <stdlib.h>
#include <iostream>

//using namespace std;

MeasuredPoints::MeasuredPoints(void)
{
}

MeasuredPoints::MeasuredPoints(std::string ad){
	file.open(ad.c_str(), std::ios::in | std::ios::binary);

	if (!file) {
		std::cerr << "Fehler beim Oeffnen Measured Points" << std::endl;
			exit(1);
	}

	this->init();		

}

MeasuredPoints::~MeasuredPoints(void)
{
}




void MeasuredPoints::init() {
	
	std::string bufferString;
	unsigned int bufferInt;

	getline(file,bufferString);
    readLevel();

	this->levelSizes.resize(maxLevel);
	this->points.resize(maxLevel);

	for (uint i=0; i<maxLevel; i++) {
		getline(file,bufferString);
		bufferInt = atoi(bufferString.c_str()); 

		this->levelSizes[i]=bufferInt;

		this->points[i].resize(levelSizes[i]);
		if(levelSizes[i] != 0) {
			for ( uint j = 0; j < levelSizes[i]; j++) {
				getline(file,bufferString);
				bufferInt = atoi(bufferString.c_str()); 
				this->points[i][j]=bufferInt;
			}
		}


	}



}