#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "VirtualFluidsBasics_EXPORT.h"

class StringUtil
{
public:
    static VirtualFluidsBasics_EXPORT std::string findAndReplace(const std::string &source, const std::string& find, const std::string& replace);
    static VirtualFluidsBasics_EXPORT std::string makeUpper(const std::string& instring);
    static VirtualFluidsBasics_EXPORT std::string makeLower(const std::string& instring);
    static VirtualFluidsBasics_EXPORT bool contains(const std::string& source, const char *find);
    static VirtualFluidsBasics_EXPORT std::string pad(const std::string& input, char pad, int length);
    static VirtualFluidsBasics_EXPORT std::string trim(const std::string &input, const std::string &trim = std::string(" \t\n"));
    static VirtualFluidsBasics_EXPORT int toInt(const std::string &input);
    static VirtualFluidsBasics_EXPORT float toFloat(const std::string &input);
    static VirtualFluidsBasics_EXPORT double toDouble(const std::string &input);
    static VirtualFluidsBasics_EXPORT bool toBool(const std::string &input);
    static VirtualFluidsBasics_EXPORT std::vector<int> toVector(const std::string& s);
    template<typename T>
    static VirtualFluidsBasics_EXPORT std::string toString(const T& t);

    static VirtualFluidsBasics_EXPORT bool endsWith(const std::string &input, const std::string &end);

private:
    StringUtil() {};
    StringUtil(const StringUtil&) {};
    virtual ~StringUtil() {};

    static bool toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base&));
};

#endif

