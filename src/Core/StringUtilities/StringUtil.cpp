#include "StringUtil.h"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

std::string StringUtil::findAndReplace(const std::string &source, const std::string& find, const std::string& replace)
{
    std::string output = source;
    size_t j;
    for (; (j = output.find(find)) != std::string::npos;)
        output.replace(j, find.length(), replace);
    return output;
}

std::string StringUtil::makeUpper(const std::string& instring)
{
    std::string output = instring;
    transform(output.begin(), output.end(), output.begin(), ::toupper);
    return output;
}

std::string StringUtil::makeLower(const std::string& instring)
{
    std::string output = instring;
    transform(output.begin(), output.end(), output.begin(), ::tolower);
    return output;
}

bool StringUtil::contains(const std::string& source, const char *find)
{
    return (0 != strstr(source.c_str(), find));
}

std::string StringUtil::pad(const std::string& input, char pad, int length)
{
    std::string outstring = input;
    for (int i = (int)outstring.length(); i < length; ++i)
        outstring += pad;
    return outstring;
}

std::string StringUtil::trim(const std::string &input, const std::string &trim /*= std::string(" \t\n")*/)
{
    if (input.size() == 0)
        return input;
    std::string temp = "";
    std::string::size_type begpos = input.find_first_not_of(trim);
    if (begpos == std::string::npos)
    {
        return temp;
    }
    else
    {
        std::string::size_type endpos = input.find_last_not_of(trim);
        temp = input.substr(begpos, endpos - begpos + 1);
    }
    return temp;
}

int StringUtil::toInt(const std::string &input)
{
    return std::stoi(input);
}

float StringUtil::toFloat(const std::string &input)
{
    return std::stof(input);
}

double StringUtil::toDouble(const std::string &input)
{
    return std::stod(input);
}

bool StringUtil::toBool(const std::string &input)
{
    bool b = 0;
    std::string trimmedInput = trim(input);
    if (!toBool(b, trimmedInput, std::boolalpha))
        throw "StringUtils::toBool() - Not a bool: " + trimmedInput;
    return b;
}

bool StringUtil::toBool(bool &t, const std::string &input, std::ios_base &(*f)(std::ios_base&))
{
    std::istringstream iss(input);
    return !(iss >> f >> t).fail();
}

std::vector<int> StringUtil::toIntVector(const std::string& input)
{
    std::vector<int> v;
    std::vector<std::string> inputEntries;
    boost::algorithm::split(inputEntries, input, boost::is_any_of("\t\n "));
    BOOST_FOREACH(std::string entry, inputEntries)
        if (entry != "")
            v.push_back(toInt(entry));
    return v;
}

std::vector<unsigned int> StringUtil::toUintVector(const std::string & input)
{
	std::vector<unsigned int> v;
	std::vector<std::string> inputEntries;
	boost::algorithm::split(inputEntries, input, boost::is_any_of("\t\n "));
	BOOST_FOREACH(std::string entry, inputEntries)
		if (entry != "")
			v.push_back(toInt(entry));
	return v;
}

std::vector<bool> StringUtil::toBoolVector(const std::string & input)
{
	std::vector<bool> v;
	std::vector<std::string> inputEntries;
	boost::algorithm::split(inputEntries, input, boost::is_any_of("\t\n "));
	BOOST_FOREACH(std::string entry, inputEntries)
	{
		bool b = 0;
		std::string trimmedInput = trim(input);
		if (toBool(b, trimmedInput, std::noboolalpha))
			v.push_back(b);
	}
	return v;
}

std::vector<std::string> StringUtil::toStringVector(const std::string & input)
{
	std::vector<std::string> v;
	std::vector<std::string> inputEntries;
	boost::algorithm::split(inputEntries, input, boost::is_any_of("\t\n "));
	BOOST_FOREACH(std::string entry, inputEntries)
		if (entry != "")
			v.push_back(toString(entry));
	return v;
}

VF_PUBLIC std::vector<double> StringUtil::toDoubleVector(const std::string & input)
{
	std::vector<double> v;
	std::vector<std::string> inputEntries;
	boost::algorithm::split(inputEntries, input, boost::is_any_of("\t\n "));
	BOOST_FOREACH(std::string entry, inputEntries)
		if (entry != "")
			v.push_back(toDouble(entry));
	return v;
}

template<typename T>
std::string StringUtil::toString(const T& t)
{
    std::ostringstream stream;
    stream << t;
    return stream.str();
}

template VF_PUBLIC std::string StringUtil::toString<int>(const int& t);


bool StringUtil::endsWith(const std::string &input, const std::string &end)
{
    return boost::algorithm::ends_with(input, end);
}