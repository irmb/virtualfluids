#ifndef BASICS_CONFIGURATIONFILE_H
#define BASICS_CONFIGURATIONFILE_H

#include <logger/Logger.h>

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <basics/utilities/UbException.h>

//! \brief  Simple configuration file
//! \details The Configuration class presented here can read and keep values of any configuration file written in a format
//! like this:
//! #
//! # Simulation parameters
//! #
//!
//! nbDimensions    = 2
//! temperature     = 25.001
//! epsilon         = 1.013e-14
//! writeLogFile    = false      # NOTE: Set to "true" in debug mode only.
//!                             #       Logging slows down the program.
//! errorMessage    = the simulation failed
//! origin          = 0.0 0.0 0.0 # x, y, z of origin
//!
//! Example how to use it:
//!
//! ConfigurationFile   config;
//! config.load(configname);
//!
//! int            nbDimensions = config.getValue<int>("nbDimensions");
//! float          temperature  = config.getValue<float>("temperature");
//! double         epsilon      = config.getValue<double>("epsilon");
//! bool           writeLogFile = config.getValue<bool>("writeLogFile");
//! string         errorMessage = config.getValue<string>("errorMessage");
//! vector<double> origin       = config.getVector<double>("origin");
//!
//! \author  Konstantin Kutscher

namespace vf::basics
{

template <class T>
T convert_to(const std::string& value)
{
    std::istringstream stream(value);
    T typedValue;
    stream >> typedValue;
    if (stream.fail())
        throw UbException(UB_EXARGS, " cannot convert \"" + value + "\" to type <" +
                                         static_cast<std::string>(typeid(typedValue).name()) + ">");

    return typedValue;
}

template <>
bool convert_to<bool>(const std::string& value);

class ConfigurationFile
{
public:
    //! clear all values
    void clear();

    //! load a configuration file
    bool load(const std::string& File);

    //! check if value associated with given key exists
    bool contains(const std::string& key) const;

    //! get vector with key
    template <class T>
    std::vector<T> getVector(const std::string& key) const;

    //! get value with key
    template <class T>
    T getValue(const std::string& key) const;

    //! get value with key and default value
    template <class T>
    T getValue(const std::string& key, T defaultValue) const;

    //! the container is public to test this class
    std::map<std::string, std::string> data;

private:
    //! the container
    //! get string with key
    std::string getValue(const std::string& key) const;

    //! remove leading and trailing tabs and spaces
    static std::string trim(const std::string& str);

    //! convert string to data type T
    // template <class T>
    // T convert_to(const std::string& value) const;

    void split(std::vector<std::string>& lst, const std::string& input, const std::string& separators,
               bool remove_empty = true) const;
};

//////////////////////////////////////////////////////////////////////////
template <class T>
std::vector<T> ConfigurationFile::getVector(const std::string& key) const
{
    std::string string_value = getValue(key);
    std::vector<T> values;
    std::vector<std::string> string_vector;
    split(string_vector, string_value, "\t\n\r;, ");
    for (std::vector<std::string>::iterator it = string_vector.begin(); it != string_vector.end(); ++it) {
        if (!(*it).empty()) {
            values.push_back(convert_to<T>(*it));
        }
    }
    return values;
}

//////////////////////////////////////////////////////////////////////////
template <class T>
T ConfigurationFile::getValue(const std::string& key) const
{
    std::string value = getValue(key);

    return convert_to<T>(value);
}

template <class T>
T ConfigurationFile::getValue(const std::string& key, T defaultValue) const
{
    if (contains(key)) {
        return getValue<T>(key);
    }
    return defaultValue;
}

ConfigurationFile loadConfig(int argc, char* argv[], std::string configPath = "./config.txt");

} // namespace vf::basics

#endif
