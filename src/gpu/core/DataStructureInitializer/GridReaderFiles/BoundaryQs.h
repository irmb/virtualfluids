#ifndef BoundaryQs_H
#define BoundaryQs_H

#include <vector>
#include <string>
#include <fstream>
#include <memory>

#include "Calculation/Calculation.h"

class Parameter;

class BoundaryQs
{
public:
    BoundaryQs();
    BoundaryQs(std::string q, bool binaer);
    BoundaryQs(std::string q, std::shared_ptr<Parameter> para, std::string str, bool binaer);
    ~BoundaryQs(void);

public:
    unsigned int getSize(unsigned int level);
    unsigned int getLevel();


private:    
    void checkFileStatus(std::string path);
    void init();
    void resizeVectors();
    void resizeVectorsPerLevel(unsigned int level, std::vector<uint32_t> &vec1D_code);

    void init_Binary();

public:
    void setIndexInVector(std::vector<std::vector<int>> &data, unsigned int level) const;
    void setValuesInVector(std::vector<std::vector<std::vector<real>>> &q27, unsigned int level) const;
    void setIndex(int *indices, unsigned int level) const;
    void setValues(real** q27, unsigned int level) const;
    void getQs(std::vector<std::vector<std::vector<real> > > &qs);
    void getIndices(std::vector<std::vector<uint> > &indices);
    //void initArray(real* ptr, unsigned int level, unsigned int column);
    //void initIndex(int *ptr, unsigned int level);

private:
    std::vector< std::vector<std::vector<real> > >values;
    std::vector< std::vector<unsigned int> >indices;

    std::ifstream file;
    unsigned int maxLevel;
    std::vector<unsigned int> levelSizes;

};

#endif