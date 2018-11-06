#ifndef Input_H
#define Input_H

#include <VirtualFluidsDefinitions.h>


#include <string>
#include <memory>
#include <istream>

namespace input
{
	class Input
	{
	public:
        static VF_PUBLIC std::unique_ptr<Input> makeInput(std::istream &stream, const std::string &inputType);

        virtual bool hasValue(const std::string &key) const = 0;
        virtual std::string getValue(const std::string &key) = 0;
	};
}

#endif
