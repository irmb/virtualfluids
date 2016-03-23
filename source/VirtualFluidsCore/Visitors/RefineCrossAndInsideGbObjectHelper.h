#ifndef RefineCrossAndInsideGbObjectHelper_H
#define RefineCrossAndInsideGbObjectHelper_H

#include <vector>

#include <Grid3D.h>
#include <GbObject3D.h>

//! \brief Refine blocks on base of bounding boxes.
//! \details You need to use <i>addGbObject()</i> to add corresponding bounding boxes. Then call <i>refine()</i>.
//! \author K. Kucher
class RefineCrossAndInsideGbObjectHelper
{
public:
   //! Constructor
   //! \param grid a smart pointer to the grid object
   //! \param maxRefineLevel an integer for maximal refinement level
   RefineCrossAndInsideGbObjectHelper(Grid3DPtr grid, int maxRefineLevel);
   virtual ~RefineCrossAndInsideGbObjectHelper(void);
   //! add geometric object
   //! \param object a smart pointer to bounding box
   //! \param refineLevel a value of refinement level for corresponding bounding box
   void addGbObject(GbObject3DPtr object, int refineLevel);
   //! start refinement
   void refine();
private:
   Grid3DPtr grid;
   std::vector<GbObject3DPtr> objects;
   std::vector<int> levels;
   int maxRefineLevel;
};

#endif 
