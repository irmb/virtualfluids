#include "MathematicaListOfListsImp.h"

#include <iomanip>
#include <limits>

std::shared_ptr<MathematicaListOfLists> MathematicaListOfListsImp::getNewInstance(std::string listName, std::vector<std::vector<double>> listOfLists)
{
    return std::shared_ptr<MathematicaListOfLists>(new MathematicaListOfListsImp(listName, listOfLists));
}

std::string MathematicaListOfListsImp::getListName()
{
    return listName;
}

MathematicaListOfListsImp::MathematicaListOfListsImp(std::string listName, std::vector<std::vector<double>> listOfLists) : listName(listName), listOfLists(listOfLists)
{
    mathematicaFunction << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    mathematicaFunction << listName << "= {";

    for (int i = 0; i < listOfLists.size(); i++) {
        mathematicaFunction << "{";
        for (int j = 0; j < listOfLists.at(i).size(); j++){
            if (j > 0)
                mathematicaFunction << ", ";
            mathematicaFunction << listOfLists.at(i).at(j);
        }
        mathematicaFunction << "}";
        if(i < listOfLists.size() - 1)
            mathematicaFunction << ", ";
    }
    mathematicaFunction << "};";
}

MathematicaListOfListsImp::MathematicaListOfListsImp()
{
}
