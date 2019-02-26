#include "VectorHelper.h"

void VectorHelper::fillVector(std::vector<int> &vector, int insertNumber) {
	fill(vector.begin(), vector.end(), insertNumber);
}

void VectorHelper::fillVector(std::vector<std::vector<int> > &vector, int insertNumber) {
	for (unsigned int i = 0; i < vector.size(); i++) {
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
	for (unsigned int i = 0; i < vector.size(); i++) {
		for (unsigned int j = 0; j < vector[i].size(); j++) {
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
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // set output default white 7;
}

void VectorHelper::dispVectorColour(const std::vector<std::vector<int>>& vector)
{
	for (unsigned int i = 0; i < vector.size(); i++) {
		for (unsigned int j = 0; j < vector[i].size(); j++) {
			makeVectorOutputColourful(vector[i][j]);
			std::cout << std::setw(4) << vector[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // set output default white 7;
}

void VectorHelper::makeVectorOutputColourful(int outputNumber)
{
	switch (outputNumber) {
	case -1:
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 8); //set output dark grey 8, dark blue 1, black 0;
		break;
	case 0:
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12); //set output bright green 10, bright red 12;
		break;
	case -5:
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 0); //set output dark grey 8, dark blue 1, black 0;
		break;
	default:
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 10); //set output bright green 10, bright red 12;
	}

}