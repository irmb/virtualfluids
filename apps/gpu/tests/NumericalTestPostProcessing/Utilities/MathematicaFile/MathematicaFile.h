#ifndef MATHEMATICA_FILE_H
#define MATHEMATICA_FILE_H

#include <fstream>
#include <memory>
#include <string>

class MathematicaFunction;

class MathematicaFile
{
public:
	static std::shared_ptr<MathematicaFile> getNewInstance(std::string filePath);

	void addMathematicaFunction(std::shared_ptr<MathematicaFunction> aMathFunc);
	void finishFile();
	

private:
	MathematicaFile();
	MathematicaFile(std::string filePath);

	std::string calcDateAndTime();

	std::string filePathTxtFile, filePathMathematicaFile;
	bool fileFinished;
	std::ofstream myFile;
};
#endif