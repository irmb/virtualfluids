#include "ConfigInput.h"
#include <errno.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>

#define COMMENT '#'

namespace input
{
    // Trim the given characters from the beginning and end of a string.
    // the default is to trim whitespace. If the string is empty or contains
    // only the trim characters, an empty string is returned.
    std::string trim(const std::string &instring, const std::string &trimstring = std::string(" \t\n"))
    {
        if (trimstring.size() == 0)
            return instring;
        std::string temp = "";
        std::string::size_type begpos = instring.find_first_not_of(trimstring);
        if (begpos == std::string::npos)
        {
            return temp;
        }
        else
        {
            std::string::size_type endpos = instring.find_last_not_of(trimstring);
            temp = instring.substr(begpos, endpos - begpos + 1);
        }
        return temp;
    }

    ConfigInput::ConfigInput(std::istream &stream) : stream(stream)
    {
        while (!stream.eof())
            this->setTokenValuePair();
    }

    ConfigInput::~ConfigInput()
    {

    }

    bool ConfigInput::hasValue(const std::string &key) const
    {
        bool valueFound = false;
        std::string keyCopy = key;
        this->makeLower(keyCopy);
        if (configEntries.find(keyCopy.c_str()) != configEntries.end())
            valueFound = true;

        return valueFound;
    }

    std::string ConfigInput::getValue(const std::string &key)
    {
        std::string keyCopy = key;
        this->makeLower(keyCopy);
        if (configEntries.find(keyCopy.c_str()) != configEntries.end())
            return (*configEntries.find(keyCopy.c_str())).second;
        return "";
    }

    //////////////////////////////////////////////////////////////////////////
    //                        private methods                               //
    //////////////////////////////////////////////////////////////////////////

    void ConfigInput::makeLower(std::string &value) const
    {
        for (unsigned i = 0; i < value.size(); i++)
            value[i] = tolower(value[i]);
    }

    void ConfigInput::setTokenValuePair()
    {
        this->eatWhiteAndComments(true);

        std::string token;
        if(!this->setToken(token))
            return;

        std::string value;
        this->setValue(value);

        configEntries.insert(String_Pair(token, value));
    }

    bool ConfigInput::setToken(std::string &token)
    {
        char tokenChar[1024];
        bool foundEqualSign = false;
        int charIndex = 0;

        this->findToken(foundEqualSign, tokenChar, charIndex);

        if (!isToken(charIndex, foundEqualSign))
            return false;

        this->nullTerminate(tokenChar, charIndex);
        token = tokenChar;
        makeLower(token);
        return true;
    }

    void ConfigInput::findToken(bool &foundEqualSign, char *token, int &i)
    {
        char ch;
        while (!(stream.get(ch)).fail())
        {
            if ((ch != '\t'))
            {
                if ((ch == '=') || (ch == ' ') || (ch == '\n') || (ch == '\r') || (ch == '\t'))
                {
                    foundEqualSign = true;
                    break;
                }
                token[i++] = ch;
            }
        }
    }


    bool ConfigInput::isToken(int charIndex, bool foundEqualSign)
    {
        if (charIndex == 0)
        {
            configEntries.insert(String_Pair("", ""));
            return false;
        }

        if (!foundEqualSign && !advanceToEqualSignOnLine())
        {
            configEntries.insert(String_Pair("", ""));
            return false;
        }
        return true;
    }

    void ConfigInput::setValue(std::string &value)
    {
        int charIndex = 0;
        char valueChar[1024];
        this->findValue(charIndex, valueChar);

        if (charIndex == 0)
            value = "";
        else
        {
            this->nullTerminate(valueChar, charIndex);
            value = valueChar;
            value = trim(value);
            this->stripLeadingAndTrailingQuotes(value);
        }
    }


    int ConfigInput::findValue(int &charIndex, char * value)
    {
        char ch;
        char c = eatWhiteAndComments(false);
        if (c != '\n')
        {
            charIndex = 0;
            while (!(stream.get(ch)).fail())
            {
                if ((ch == '\t') || (ch == '\r') || (ch == '\n') || (ch == '#'))
                {
                    while (ch != '\n')
                    {
                        if (stream.get(ch).fail()) break;
                    }
                    break;
                }
                else
                {
                    value[charIndex++] = ch;
                }
            }
        }       
        return charIndex;
    }



    void ConfigInput::stripLeadingAndTrailingQuotes(std::string &m_value)
    {
        if (m_value[0] == '"')
            m_value = m_value.substr(1);
        if (m_value[m_value.length() - 1] == '"')
            m_value = m_value.substr(0, m_value.length() - 1);
    }

    void ConfigInput::nullTerminate(char *value, int &i)
    {
        value[i++] = '\0';
    }

    bool ConfigInput::advanceToEqualSignOnLine()
    {
        char ch;
        bool foundEqual = false;
        while (!(stream.get(ch)).fail())
        {
            if (isNewLine(ch) || isCarriageReturn(ch)) 
                break;
            if (isEqualSign(ch))
            {
                foundEqual = true;
                break;
            }
        }
        return foundEqual;
    }

    char ConfigInput::eatWhiteAndComments(bool traverseNewlines)
    {
        char ch;
        bool isComment = false;

        while (!(stream.get(ch)).fail())
            if (isCommentSign(ch))
                isComment = true;
            else if (isNewLine(ch))
            {
                isComment = false;
                if (!traverseNewlines)
                    return(ch);
            }
            else if (isRegularChar(isComment, ch))
            {
                stream.putback(ch);
                return 0;
            }
            return 0;
    }

    bool ConfigInput::isRegularChar(bool isComment, char ch)
    {
        return (!isComment) && (ch != ' ') && (ch != '\t') && (ch != '\r');
    }

    bool ConfigInput::isCommentSign(char ch)
    {
        return ch == COMMENT;
    }

    bool ConfigInput::isEqualSign(char ch)
    {
        return ch == '=';
    }

    bool ConfigInput::isNewLine(char ch)
    {
        return ch == '\n';
    }

    bool ConfigInput::isCarriageReturn(char ch)
    {
        return ch == '\r';
    }


}
