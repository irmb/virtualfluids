//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEINPUTASCII_H
#define UBFILEINPUTASCII_H

#include <fstream>
#include <iostream>
#include <string>

#include "basics_export.h"

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>

/*=========================================================================*/
/*  UbFileInputASCII                                                       */
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

class BASICS_EXPORT UbFileInputASCII : public UbFileInput
{
public:
    UbFileInputASCII() : UbFileInput() {}
    UbFileInputASCII(std::string filename);

    bool open(std::string filename) override;

    std::string getFileName() override;
    void skipLine() override; // Springt zur naechsten Zeile

    void readLine() override;
    std::string readStringLine() override;
    int readInteger() override; // Liest einen Int-Wert ein
    long long readLongLong();   // Liest einen long-Wert ein

    std::size_t readSize_t() override;
    double readDouble() override;                 // Liest einen double-Wert ein
    float readFloat() override;                   // Liest einen float-Wert ein
    bool readBool() override;                     // Liest einen bool-Wert ein
    char readChar() override;                     // Liest einen char-Wert ein
    std::string readString() override;            // Liest ein Wort ein
    std::string readLineTill(char stop) override; // Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
    std::string parseString() override;

    bool containsString(const std::string &var) override;
    void setPosAfterLineWithString(const std::string &var) override;
    int readIntegerAfterString(const std::string &var) override;
    double readDoubleAfterString(const std::string &var) override;
    bool readBoolAfterString(const std::string &var) override;
    std::string readStringAfterString(const std::string &var) override;

    FILETYPE getFileType() override { return ASCII; }

    template <typename T>
    friend inline UbFileInputASCII &operator>>(UbFileInputASCII &file, T &data)
    {
        file.infile >> data;
        return file;
    }
};

#endif // UBFILEINPUTASCII_H
