#include "MathematicaFile.h"

#include "Utilities/MathematicaFunction/MathematicaFunktion.h"

#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>

std::shared_ptr<MathematicaFile> MathematicaFile::getNewInstance(std::string filePath)
{
    return std::shared_ptr<MathematicaFile>(new MathematicaFile(filePath));
}

void MathematicaFile::addMathematicaFunction(std::shared_ptr<MathematicaFunction> aMathFunc)
{
    if (!fileFinished)
        myFile << aMathFunc->getFunction() << std::endl;
}

void MathematicaFile::finishFile()
{
    if (!fileFinished) {
        fileFinished = true;
        myFile.close();

        std::filesystem::rename(filePathTxtFile, filePathMathematicaFile);
        std::system(filePathMathematicaFile.c_str());
    }
}

MathematicaFile::MathematicaFile(std::string filePath)
{
    fileFinished = false;

    std::ostringstream basicFilePath, textFile, mathematicaFile;
    basicFilePath << filePath << "/NumericalTestPostProcessing_" << calcDateAndTime();
    textFile << basicFilePath.str() << ".txt";
    mathematicaFile << basicFilePath.str() << ".nb";
    filePathTxtFile = textFile.str();
    filePathMathematicaFile = mathematicaFile.str();

    myFile.open(filePathTxtFile);
}

std::string MathematicaFile::calcDateAndTime()
{
    std::ostringstream oss;
    time_t now = time(NULL);
    struct tm nowLocal = *localtime(&now);
    oss << std::setfill('0') << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
    return oss.str();
}

MathematicaFile::MathematicaFile()
{
}
