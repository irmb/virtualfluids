//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEINPUTBINARY_H
#define UBFILEINPUTBINARY_H

#include <fstream>
#include <iostream>
#include <string>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbFileInputBinary                                                      */
/*                                                                         */
/**
...
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 23.11.04
*/

/*
usage: ...
*/

class UbFileInputBinary : public UbFileInput
{
public:
    UbFileInputBinary() : UbFileInput() {}
    UbFileInputBinary(std::string filename);

    bool open(std::string filename) override;

    void skipLine() override; // Springt zur naechsten Zeile
    void readLine() override;
    std::string readStringLine() override;
    std::size_t readSize_t() override;
    int readInteger() override;                   // Liest einen Int-Wert ein
    double readDouble() override;                 // Liest einen double-Wert ein
    float readFloat() override;                   // Liest einen float-Wert ein
    bool readBool() override;                     // Liest einen bool-Wert ein
    char readChar() override;                     // Liest einen char-Wert ein
    std::string readString() override;            // Liest ein Wort ein
    std::string readLineTill(char stop) override; // Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
    std::string parseString() override;           // Liest

    bool containsString(const std::string &var) override;
    void setPosAfterLineWithString(const std::string &var) override;
    int readIntegerAfterString(const std::string &var) override;
    double readDoubleAfterString(const std::string &var) override;
    bool readBoolAfterString(const std::string &var) override;
    std::string readStringAfterString(const std::string &var) override;

    FILETYPE getFileType() override { return BINARY; }

    template <typename T>
    friend inline UbFileInputBinary &operator>>(UbFileInputBinary &file, T &data)
    {
        file.infile.read((char *)&data, sizeof(T));
        return file;
    }

    template <typename T>
    void readVector(std::vector<T> &v)
    {
        size_t size = v.size();
        if (size > 0) {
            infile.read((char *)&v[0], sizeof(T) * size);
        }
    }
};

#endif
