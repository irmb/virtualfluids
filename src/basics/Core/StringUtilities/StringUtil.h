#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "basics_export.h"

#define SSTR(x) static_cast<std::ostringstream &>((std::ostringstream() << std::dec << x)).str()

class StringUtil
{
public:
    static BASICS_EXPORT std::string findAndReplace(const std::string &source, const std::string &find,
                                                    const std::string &replace);
    static BASICS_EXPORT std::string makeUpper(const std::string &instring);
    static BASICS_EXPORT std::string makeLower(const std::string &instring);
    static BASICS_EXPORT std::vector<std::string> split(const std::string &input, const std::string &delim = " ");
    static BASICS_EXPORT bool contains(const std::string &source, const char *find);
    static BASICS_EXPORT std::string pad(const std::string &input, char pad, int length);
    static BASICS_EXPORT std::string trim(const std::string &input, const std::string &trim = std::string(" \t\n"));
    static BASICS_EXPORT int toInt(const std::string &input);
    static BASICS_EXPORT float toFloat(const std::string &input);
    static BASICS_EXPORT double toDouble(const std::string &input);
    static BASICS_EXPORT bool toBool(const std::string &input);
    static BASICS_EXPORT std::vector<int> toIntVector(const std::string &s);
    static BASICS_EXPORT std::vector<unsigned int> toUintVector(const std::string &s);
    static BASICS_EXPORT std::vector<bool> toBoolVector(const std::string &s);
    static BASICS_EXPORT std::vector<std::string> toStringVector(const std::string &s);
    static BASICS_EXPORT std::vector<double> toDoubleVector(const std::string &s);
    template <typename T>
    static BASICS_EXPORT std::string toString(const T &t);

    static BASICS_EXPORT bool endsWith(const std::string &input, const std::string &end);

private:
    StringUtil() = default;
    ;
    StringUtil(const StringUtil &) = default;
    ;
    virtual ~StringUtil() = default;
    ;

    static bool toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base &));
};

#endif
