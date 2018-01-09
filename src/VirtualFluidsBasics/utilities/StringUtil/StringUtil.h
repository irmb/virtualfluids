#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include <VirtualFluidsDefinitions.h>

class StringUtil
{
public:
    static VF_PUBLIC std::string findAndReplace(const std::string &source, const std::string& find, const std::string& replace);
    static VF_PUBLIC std::string makeUpper(const std::string& instring);
    static VF_PUBLIC std::string makeLower(const std::string& instring);
    static VF_PUBLIC bool contains(const std::string& source, const char *find);
    static VF_PUBLIC std::string pad(const std::string& input, char pad, int length);
    static VF_PUBLIC std::string trim(const std::string &input, const std::string &trim = std::string(" \t\n"));
    static VF_PUBLIC int toInt(const std::string &input);
    static VF_PUBLIC float toFloat(const std::string &input);
    static VF_PUBLIC double toDouble(const std::string &input);
    static VF_PUBLIC bool toBool(const std::string &input);
    static VF_PUBLIC std::vector<int> toVector(const std::string& s);
    template<typename T>
    static VF_PUBLIC std::string toString(const T& t);

    static VF_PUBLIC bool endsWith(const std::string &input, const std::string &end);

private:
    StringUtil() {};
    StringUtil(const StringUtil&) {};
    virtual ~StringUtil() {};

    static bool toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base&));
};

#endif

