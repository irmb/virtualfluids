#ifndef RefineAroundGbObjectHelper_H
#define RefineAroundGbObjectHelper_H

#include <vector>

#include <Grid3D.h>
#include <D3Q27TriFaceMeshInteractor.h>

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
   RefineAroundGbObjectHelper(Grid3DPtr grid, int maxRefineLevel, D3Q27TriFaceMeshInteractorPtr objectIter, double startDistance, double stopDistance, CommunicatorPtr comm);
   virtual ~RefineAroundGbObjectHelper(void);
   //! start refinement
   void refine();
private:
   Grid3DPtr grid;
   D3Q27TriFaceMeshInteractorPtr objectIter;
   int refineLevel;
   double startDistance, stopDistance;
   CommunicatorPtr comm;
};

#endif 
