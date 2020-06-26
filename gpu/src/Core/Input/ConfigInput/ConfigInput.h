#ifndef ConfigInput_H
#define ConfigInput_H
#include <string>
#include <vector>
#include <istream>
#include <memory>
#include <map>
#include <list>

#include <VirtualFluidsDefinitions.h>
#include "../Input.h"

namespace input 
{
    class  ConfigInput : public Input
    {
    public:
        VF_PUBLIC ConfigInput(std::istream &stream);
        VF_PUBLIC virtual ~ConfigInput(void);
   
        VF_PUBLIC bool hasValue(const std::string &key) const;
        VF_PUBLIC std::string getValue(const std::string &key);

    protected:
        virtual void setTokenValuePair();
        void setValue(std::string &value);
        bool setToken(std::string &token);
        bool isToken(int charIndex, bool foundEqualSign);
        int findValue(int &charIndex, char * value);
        void stripLeadingAndTrailingQuotes(std::string &m_value);
        void nullTerminate(char * token, int &i);
        void findToken(bool &foundEqualSign, char * token, int &i);
        char eatWhiteAndComments(bool traverse_newlines = true);
        bool isRegularChar(bool isComment, char ch);
        bool isNewLine(char ch);
        bool isCommentSign(char ch);
        bool advanceToEqualSignOnLine();
        bool isCarriageReturn(char ch);
        bool isEqualSign(char ch);
        void makeLower(std::string &instring) const;

    protected:
        std::istream &stream;
        typedef std::pair <std::string, std::string> String_Pair;
        std::map<std::string, std::string> configEntries;
    };
}
#endif