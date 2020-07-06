#include "ConfigFile.h"
#include "Utilities/StringUtil.hpp"
#include <errno.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>


ConfigFile::ConfigFile( const char *strConfigFile ) 
{
   if ( strConfigFile )
      m_sConfigFile = strConfigFile;
}
ConfigFile::~ConfigFile()
{
}
bool ConfigFile::read() 
{
   m_in.open(m_sConfigFile.c_str(),std::ios::in);
   if (m_in.fail())
   {
      return false;
   }
   while (!m_in.eof())
   {
      //--------------------------------------------------------
      // Get a token and value.
      // This gives values to member vars: m_token and m_value.
      //----------------------------------------------------------
      get_token_and_value();
      if ( m_token.length() )
         m_ConfigEntries.insert( String_Pair(m_token, m_value) );
   }
   m_in.close();
   return true;
}
void ConfigFile::get_token_and_value(void)
{
   char token[1024];
   char ch;
   bool found_equal=false;
   int i=0;
   eat_white_and_comments();
   while(!(m_in.get(ch)).fail())
   {
      if ((ch != '\t'))
      {
         if ( (ch == '=') || (ch == ' ') || (ch == '\n') || (ch == '\r') || 
            (ch == '\t'))
         {
            if (ch == '=')found_equal=true;
            break;
         }
         token[i++]=ch;
      }
   }
   if (i==0)
   {
      // It didn’t find a token, in this case.
      m_token="";
      m_value="";
      return;
   }
   // Null-terminate the token that was found.
   token[i++]='\0';
   m_token = token;
   makeLower(m_token);
   // Advance to the equal sign, if need be.
   if (!found_equal)
   {
      if (!advance_to_equal_sign_on_line())
      {
         // The token had no value.
         m_token="";
         m_value="";
         return;
      }
   }
   // Get the token’s value.
   i=0;
   char c = eat_white_and_comments(false);
   if ( c != '\n' )
   {
      i=0;
      while(!(m_in.get(ch)).fail())
      {
         if ((ch == '\t') || (ch == '\r') ||  (ch == '\n') || (ch == '#') )
         {
            while (ch!='\n')
            {
               if (m_in.get(ch).fail()) break;
            }
            break;
         }
         else
         {
            token[i++]=ch;
         }
      }
   }
   if (i==0)
   {
      // This token had no value.
      m_value="";
   }
   else
   {
      token[i++]='\0';
      m_value=token;
      // Remove leading/trailing spaces.
      m_value = StringUtil::trim(m_value);
      // Strip leading and trailing quotes, if there are any.
      if ( m_value[0] == '"' )
         m_value = m_value.substr( 1 );
      if ( m_value[ m_value.length() -1 ] == '"' )
         m_value = m_value.substr( 0, m_value.length()-1 );
   }
}
bool ConfigFile::advance_to_equal_sign_on_line()
{
   char ch;
   bool found_equal=false;
   while ( !(m_in.get(ch)).fail() )
   {
      if ((ch=='\r')||(ch=='\n')) break;
      if (ch == '=')
      {
         found_equal=true;
         break;
      }
   }
   return found_equal;
}
char ConfigFile::eat_white_and_comments(bool traverse_newlines)
{
   char ch;
   bool in_comment;
   in_comment = false;
   while (!(m_in.get(ch)).fail())
      if (ch == '#')
         in_comment = true;
      else if (ch == '\n')
      {
         in_comment = false;
         if (!traverse_newlines)
         {
            return(ch); // Stop eating.
         }
      }
      else if ((!in_comment) && (ch != ' ') &&
         (ch != '\t') && (ch != '\r'))
      {
         m_in.putback(ch);
         return 0;
      }
      return 0;
}
void ConfigFile::makeLower(std::string &instring)
{
   for(unsigned i=0; i < instring.size();
      i++)
   {
      instring[i] = tolower(instring[i]);
   }
}
bool ConfigFile::hasValue( const char *key ) 
{
   bool bRet = false;
   std::string sKey = key;
   makeLower( sKey );
   if ( m_ConfigEntries.find( sKey.c_str() ) != m_ConfigEntries.end() )
   {
      bRet = true;
   }
   return bRet;
}
std::string ConfigFile::getValue( const char *key )
{
   std::string sKey = key;
   makeLower( sKey );
   if ( m_ConfigEntries.find( sKey.c_str() ) != m_ConfigEntries.end() )
   {
      std::map<std::string, std::string>::iterator iter;
      iter =  m_ConfigEntries.find(sKey.c_str());
      return (*iter).second;
   }
   return "";
}
void ConfigFile::setValue( const char *key, const char *value )
{
   std::string sKey = key;
   makeLower( sKey );
   m_ConfigEntries[sKey] = value;
}