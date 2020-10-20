#ifndef RefineAroundGbObjectHelper_H
#define RefineAroundGbObjectHelper_H

#include <PointerDefinitions.h>

class Grid3D;
class Communicator;
class D3Q27TriFaceMeshInteractor;

//! \brief Refine blocks on base of bounding boxes.
//! \details You need to use <i>addGbObject()</i> to add corresponding bounding boxes. Then call <i>refine()</i>.
//! \author K. Kucher
class RefineAroundGbObjectHelper
{
public:
    //! Constructor
    //! \param grid a smart pointer to the grid object
    //! \param maxRefineLevel an integer for maximal refinement level
    //! \param objectIter a D3Q27TriFaceMeshInteractor object - represent geometry which should be refinement
    //! \param startDistance start distance from geometry for refinement
    //! \param stopDistance stop distance from geometry for refinement
    RefineAroundGbObjectHelper(SPtr<Grid3D> grid, int maxRefineLevel, SPtr<D3Q27TriFaceMeshInteractor> objectIter,
                               double startDistance, double stopDistance, SPtr<Communicator> comm);
    virtual ~RefineAroundGbObjectHelper();
    //! start refinement
    void refine();

private:
    SPtr<Grid3D> grid;
    SPtr<D3Q27TriFaceMeshInteractor> objectIter;
    int refineLevel;
    double startDistance, stopDistance;
    SPtr<Communicator> comm;
};

#endif
