#ifndef CONFIGFILE_H
#define CONFIGFILE_H
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <list>

class ConfigFileReader
{
public:
   ConfigFileReader(const char *strFileName);
   virtual ~ConfigFileReader(void);
   bool read(void);
   bool hasValue( const char *key );
   std::string getValue( const char *key );
   void setValue( const char *key, const char *value );
protected:
   virtual void get_token_and_value();
   virtual char eat_white_and_comments(bool traverse_newlines=true);
   virtual bool advance_to_equal_sign_on_line();
   virtual void makeLower(std::string &instring);
protected:
   std::fstream m_in;
   std::string m_token;
   std::string m_value;
   std::string m_sConfigFile;
   typedef std::pair <std::string, std::string> String_Pair;
   std::map<std::string, std::string> m_ConfigEntries;
};
#endif
