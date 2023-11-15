#ifndef MATHEMATICA_LIST_OF_LISTS_IMP_H
#define MATHEMATICA_LIST_OF_LISTS_IMP_H

#include "MathematicaListOfLists.h"

#include <memory>
#include <vector>

class MathematicaListOfListsImp : public MathematicaListOfLists
{
public:

    static std::shared_ptr<MathematicaListOfLists> getNewInstance(std::string listName, std::vector<std::vector<double> > listOfLists);

    std::string getListName();

    
private:
    MathematicaListOfListsImp();
    MathematicaListOfListsImp(std::string listName, std::vector<std::vector<double> > listOfLists);

    std::string listName;
    std::vector<std::vector<double> > listOfLists;
};
#endif