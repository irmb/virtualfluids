#ifndef LOGFILE_READER_H
#define LOGFILE_READER_H

#include <memory>
#include <string>
#include <vector>

class LogFileData;

class LogFileReader
{
public:
    static std::shared_ptr<LogFileReader> getInstance();
    
    std::shared_ptr<LogFileData> readLogFileToLogFileData(std::string filePath);
    std::vector<std::shared_ptr<LogFileData> > readLogFilesInDirectoryToLogFileData(std::string directory);
    

private:
    LogFileReader();

    std::vector<std::string> getAllFilesInDir(const std::string &dirPath, const std::string &fileExtension);
    std::string removeCharsFromString(std::string str, char* charsToRemove);
    bool checkEqualDouble(double one, double two);

};

#endif