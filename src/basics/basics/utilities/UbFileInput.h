//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBFILEINPUT_H
#define UBFILEINPUT_H

#include <fstream>
#include <iostream>
#include <string>

#include <cstdlib> //atoi
#include <cstring> //strstr

#include <basics/utilities/UbException.h>

/*=========================================================================*/
/*  UbFileInput                                                            */
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

class UbFileInput
{
public:
    enum FILETYPE { ASCII, BINARY };

public:
    UbFileInput() : filename("") {}
    virtual ~UbFileInput() = default;

    virtual bool operator!() { return !(infile); }
    virtual bool isOpen() { return !(!(infile)); }

    virtual bool open(std::string filename) = 0;
    virtual void close() { infile.close(); }
    virtual int eof() { return infile.eof(); }

    virtual void skipLine()                     = 0; // Springt zur naechsten Zeile
    virtual void readLine()                     = 0;
    virtual std::string readStringLine()        = 0;
    virtual int readInteger()                   = 0; // Liest einen Int-Wert ein
    virtual std::size_t readSize_t()            = 0;
    virtual double readDouble()                 = 0; // Liest einen double-Wert ein
    virtual float readFloat()                   = 0; // Liest einen float-Wert ein
    virtual bool readBool()                     = 0; // Liest einen bool-Wert ein
    virtual char readChar()                     = 0; // Liest einen char-Wert ein
    virtual std::string readString()            = 0; // Liest ein Wort ein
    virtual std::string readLineTill(char stop) = 0; // Liest gesamte Zeile ein bis zu einem bestimmten Zeichen
    virtual std::string parseString()           = 0; // Liest

    virtual void setCommentIndicator(char commentindicator) { this->commentindicator = commentindicator; }

    virtual bool containsString(const std::string &var)               = 0;
    virtual void setPosAfterLineWithString(const std::string &var)    = 0;
    virtual int readIntegerAfterString(const std::string &var)        = 0;
    virtual double readDoubleAfterString(const std::string &var)      = 0;
    virtual bool readBoolAfterString(const std::string &var)          = 0;
    virtual std::string readStringAfterString(const std::string &var) = 0;

    virtual std::string getFileName() { return this->filename; }

    // returns file extension:
    // e.g. "./../test/ich.inp" -> "inp", "./../test/ich" -> ""
    virtual std::string getFileExtension()
    {
        std::size_t pos1 = filename.rfind("/");
        if (pos1 == std::string::npos)
            pos1 = 0;
        std::size_t pos2 = filename.rfind(".");
        if (pos2 != std::string::npos && pos2 > pos1)
            return filename.substr(pos2 + 1);

        return "";
    }

    // returns "ASCII", "BINARY"
    virtual FILETYPE getFileType() = 0;

protected:
    std::ifstream infile;
    std::string filename;
    char commentindicator{ 'C' };
};

#endif // UBFILEINPUT_H
