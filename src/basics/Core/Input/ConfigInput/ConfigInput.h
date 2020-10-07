#ifndef ConfigInput_H
#define ConfigInput_H
#include <string>
#include <vector>
#include <istream>
#include <memory>
#include <map>
#include <list>

#include "basics_export.h"
#include "../Input.h"

namespace input 
{
    class  ConfigInput : public Input
    {
    public:
        BASICS_EXPORT ConfigInput(std::istream &stream);
        BASICS_EXPORT ~ConfigInput() override;
   
        BASICS_EXPORT bool hasValue(const std::string &key) const override;
        BASICS_EXPORT std::string getValue(const std::string &key) override;

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
        using String_Pair = std::pair <std::string, std::string>;
        std::map<std::string, std::string> configEntries;
    };
}
#endif
