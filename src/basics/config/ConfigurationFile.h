#ifndef BASICS_CONFIGURATIONFILE_H
#define BASICS_CONFIGURATIONFILE_H

#include <logger/Logger.h>

#include <filesystem>
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include <basics/utilities/UbException.h>

//! \brief  Simple configuration file
//! \details The Configuration class presented here can read and keep values of any configuration file written in a format like this:
//!#
//!# Simulation parameters
//!#
//!
//!nbDimensions    = 2
//!temperature     = 25.001
//!epsilon         = 1.013e-14
//!writeLogFile    = false      # NOTE: Set to "true" in debug mode only.
//!                             #       Logging slows down the program.
//!errorMessage    = the simulation failed
//!origin          = 0.0 0.0 0.0 # x, y, z of origin
//!
//!Example how to use it:
//!
//!ConfigurationFile   config;
//!config.load(configname);
//!
//!int            nbDimensions = config.getValue<int>("nbDimensions");
//!float          temperature  = config.getValue<float>("temperature");
//!double         epsilon      = config.getValue<double>("epsilon");
//!bool           writeLogFile = config.getValue<bool>("writeLogFile");
//!string         errorMessage = config.getValue<string>("errorMessage");
//!vector<double> origin       = config.getVector<double>("origin");
//!
//! \author  Konstantin Kutscher


namespace vf::basics
{


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
   template<class T>
   std::vector<T> getVector(const std::string& key) const;

   //! get value with key
   template<class T>
   T getValue(const std::string& key) const;

   //! get value with key and default value
   template<class T>
   T getValue(const std::string& key, T defaultValue) const;

   //! the container is public to test this class
   std::map<std::string, std::string> data;

private:
   //! the container
   //! get string with key
   std::string  getValue(const std::string& key) const;

   //! remove leading and trailing tabs and spaces
   static std::string trim(const std::string& str);

   //! convert string to data type T
   template<class T>
   T convert_to(const std::string& str) const;

   void split(std::vector<std::string>& lst, const std::string& input, const std::string& separators, bool remove_empty = true) const;
};


//////////////////////////////////////////////////////////////////////////
template<class T>
std::vector<T> ConfigurationFile::getVector(const std::string& key) const
{
   std::string string_value = getValue(key);
   std::vector<T> values;
   std::vector<std::string> string_vector;
   split(string_vector, string_value, "\t\n\r;, ");
   for (std::vector<std::string>::iterator it = string_vector.begin(); it != string_vector.end(); ++it)
   {
      if (*it != "")
      {
         values.push_back(convert_to<T>(*it));
      }
   }
   return values;
}
//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::convert_to(const std::string& value) const
{
   if constexpr (std::is_same_v<T, bool>)
   {
      return (value == "true");
   }

   std::istringstream stream(value);
   T t;
   stream >> t;
    if (stream.fail())
      throw UbException(UB_EXARGS, " cannot convert \"" + value + "\" to type <" + static_cast<std::string>(typeid(t).name()) + ">");

   return t;
}

//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::getValue(const std::string& key) const
{
   std::string value = getValue(key);

   return convert_to<T>(value);
}

template<class T>
T ConfigurationFile::getValue(const std::string& key, T defaultValue) const
{
   if (contains(key))
   {
      return getValue<T>(key);
   }
   else
   {
      return defaultValue;
   }
}

static ConfigurationFile loadConfig(int argc, char *argv[], std::string configPath = "./config.txt")
{
   // the config file's default path can be replaced by passing a command line argument

   if (argc > 1)
   {
      configPath = argv[1];
      VF_LOG_INFO("Using command line argument for config path: {}", configPath);
   } else {
      VF_LOG_INFO("Using default config path: {}", configPath);
   }

   vf::basics::ConfigurationFile config;
   config.load(configPath);
   return config;
}


}

#endif
