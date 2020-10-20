//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEOUTPUTASCII_H
#define UBFILEOUTPUTASCII_H

#include <fstream>
#include <iomanip>
#include <iostream>

#include "basics_export.h"

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileOutput.h>

/*=========================================================================*/
/*  UbFileOutputASCII                                                             */
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

class BASICS_EXPORT UbFileOutputASCII : public UbFileOutput
{
public:
    UbFileOutputASCII() : UbFileOutput() {}
    UbFileOutputASCII(const std::string &filename, const bool &createPath = true, const int &precision = 15);
    UbFileOutputASCII(const std::string &filename, CREATEOPTION opt, const bool &createPath = true,
                      const int &precision = 15);

    bool open(const std::string &filename, CREATEOPTION opt = OUTFILE) override;

    void writeBool(const bool &value, const int &width = 0) override;
    void writeDouble(const double &value, const int &width = 0) override;
    void writeFloat(const float &value, const int &width = 0) override;
    void writeInteger(const int &value, const int &width = 0) override;
    void writeSize_t(const std::size_t &value, const int &width = 0) override;
    void writeChar(const char &value, const int &width = 0) override;
    void writeString(const std::string &value, const int &width = 0) override;
    void writeStringOnly(const std::string &value) override;
    void writeLine(const std::string &value, const int &width = 0) override;
    void writeLine() override;

    void setPrecision(const int &precision) override;
    int getPrecision() override { return (int)outfile.precision(); }

    void setCommentIndicator(char commentindicator) override { this->commentindicator = commentindicator; }

    void writeCommentLine(const std::string &line) override;
    void writeCommentLine(char indicator, const std::string &line) override;
    void writeCopyOfFile(const std::string &filename) override;

    FILETYPE getFileType() override { return ASCII; }

    template <typename T>
    friend inline UbFileOutputASCII &operator<<(UbFileOutputASCII &file, const T &data)
    {
        file.outfile << data;
        return file;
    }
};

#endif
