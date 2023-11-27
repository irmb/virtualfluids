#ifndef RefineCrossAndInsideGbObjectHelper_H
#define RefineCrossAndInsideGbObjectHelper_H

#include <PointerDefinitions.h>
#include <vector>

namespace vf::parallel {class Communicator;}
class Grid3D;
class GbObject3D;

//! \brief Refine blocks on base of bounding boxes.
//! \details You need to use <i>addGbObject()</i> to add corresponding bounding boxes. Then call <i>refine()</i>.
//! \author K. Kutscher
class RefineCrossAndInsideGbObjectHelper
{
public:
    //! Constructor
    //! \param grid a smart pointer to the grid object
    //! \param maxRefineLevel an integer for maximal refinement level
    RefineCrossAndInsideGbObjectHelper(SPtr<Grid3D> grid, int maxRefineLevel, std::shared_ptr<vf::parallel::Communicator> comm);
    virtual ~RefineCrossAndInsideGbObjectHelper();
    //! add geometric object
    //! \param object a smart pointer to bounding box
    //! \param refineLevel a value of refinement level for corresponding bounding box
    void addGbObject(SPtr<GbObject3D> object, int refineLevel);
    //! start refinement
    void refine();

private:
    SPtr<Grid3D> grid;
    std::vector<SPtr<GbObject3D>> objects;
    std::vector<int> levels;
    int maxRefineLevel;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif
