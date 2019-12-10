#ifndef ENUM_MAPPER_IMP_H
#define ENUM_MAPPER_IMP_H

#include <map>
#include <string>

template <typename T>
class EnumMapperImp
{
public:
	std::string getString(T enumeration);
	T getEnum(std::string name);
	void addEnum(T enumeration, std::string name);

private:
	std::map< T, std::string> enumMap;
	std::map< std::string, T> stringMap;

};
#endif 

template<typename T>
inline std::string EnumMapperImp<T>::getString(T enumeration)
{
	std::map< T, std::string>::iterator it;
	it = enumMap.find(enumeration);
	if (it == enumMap.end()) {
		throw std::exception("Enumeration is not registered.");
	}
	else
		return it->second;
}

template<typename T>
inline T EnumMapperImp<T>::getEnum(std::string name)
{
	std::map< std::string, T>::iterator it;
	it = stringMap.find(name);
	if (it == stringMap.end()) {
		throw std::exception("String is not registered.");
	}
	else
		return it->second;
}

template<typename T>
inline void EnumMapperImp<T>::addEnum(T enumeration, std::string name)
{
	enumMap.insert(std::pair< T, std::string>(enumeration, name));
	stringMap.insert(std::pair< std::string, T>(name, enumeration));
}
