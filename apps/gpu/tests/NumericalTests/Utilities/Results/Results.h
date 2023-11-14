#ifndef RESULTS_H
#define RESULTS_H

#include <vector>

class Results
{
public:
    virtual ~Results() = default;
    virtual int getNumberOfTimeSteps() = 0;
    virtual std::vector<std::vector<double> > getVx() = 0;
    virtual std::vector<std::vector<double> > getVy() = 0;
    virtual std::vector<std::vector<double> > getVz() = 0;
    virtual int getNumberOfXNodes() = 0;
    virtual int getNumberOfYNodes() = 0;
    virtual int getNumberOfZNodes() = 0;
    virtual std::vector<std::vector<double> > getXNodes() = 0;
    virtual std::vector<std::vector<double> > getYNodes() = 0;
    virtual std::vector<std::vector<double> > getZNodes() = 0;
    virtual int getTimeStepLength() = 0;
    virtual std::vector<unsigned int> getTimeSteps() = 0;
    virtual std::vector < std::vector<unsigned int> > getLevels() = 0;

    virtual bool checkYourData() = 0;

    virtual int getL0() = 0;

};
#endif