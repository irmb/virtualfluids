#include "VectorHelper.h"

#include <iostream>
#include <iomanip>	//formatting output streams

#include "ConsoleColor.h"

void VectorHelper::fillVector(std::vector<int> &vector, int insertNumber) {
	fill(vector.begin(), vector.end(), insertNumber);
}

void VectorHelper::fillVector(std::vector<std::vector<int> > &vector, int insertNumber) {
	for (uint i = 0; i < vector.size(); i++) {
		fill(vector[i].begin(), vector[i].end(), insertNumber);
	}
}

void VectorHelper::dispVector(const std::vector<int> &vector)
{
	for (int number : vector) {
		std::cout << std::setw(4) << number;
	}
	std::cout << std::endl;
}

void VectorHelper::dispVector(const std::vector<std::vector<int> > &vector)
{
	for (uint i = 0; i < vector.size(); i++) {
		for (uint j = 0; j < vector[i].size(); j++) {
			std::cout << std::setw(4) << vector[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void VectorHelper::dispVectorColour(const std::vector<int> &vector)
{
	for (int number : vector) {
		makeVectorOutputColourful(number);
		std::cout << std::setw(4) << number;
	}
	std::cout << std::endl;
	ConsoleColor::setDefaultWhite();
}

void VectorHelper::dispVectorColour(const std::vector<std::vector<int>>& vector)
{
	for (uint i = 0; i < vector.size(); i++) {
		for (uint j = 0; j < vector[i].size(); j++) {
			makeVectorOutputColourful(vector[i][j]);
			std::cout << std::setw(4) << vector[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	ConsoleColor::setDefaultWhite();
}

void VectorHelper::makeVectorOutputColourful(int outputNumber)
{
	switch (outputNumber) {
	case -1:
		ConsoleColor::setDarkGrey();
		break;
	case 0:
		ConsoleColor::setBrightRed();
		break;
	case -5:
		ConsoleColor::setBlack();
		break;
	default:
		ConsoleColor::setBrightGreen();
	}

}