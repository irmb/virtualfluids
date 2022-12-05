//! \file InterpolationCellGrouper.h
//! \ingroup GPU
//! \author Anna Wellmann
//! \details See [master thesis of Anna Wellmann]

#ifndef InterpolationCellGrouper_H
#define InterpolationCellGrouper_H

#include <memory>
#include <basics/Core/DataTypes.h>

class Parameter;
class GridBuilder;

class InterpolationCellGrouper
{
public:
    //! \brief Construct InterpolationCellGrouper object
    InterpolationCellGrouper(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder);

    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief Split the interpolation cells from coarse to fine into border an bulk
    //! \details For communication hiding, the interpolation cells from the coarse to the fine grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! see [master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")]
    void splitCoarseToFineIntoBorderAndBulk(uint level);

    //! \brief Split the interpolation cells from fine to coarse into border an bulk
    //! \details For communication hiding, the interpolation cells from the fine to the coarse grid need to be split
    //! into two groups:
    //!
    //! - cells which are at the border between two gpus --> "border"
    //!
    //! - the other cells which are not directly related to the communication between the two gpus --> "bulk"
    //!
    //! See [master thesis of Anna Wellmann (p. 62-68: "Überdeckung der reduzierten Kommunikation")]
    void splitFineToCoarseIntoBorderAndBulk(uint level);

protected:
    //////////////////////////////////////////////////////////////////////////
    // split interpolation cells
    //////////////////////////////////////////////////////////////////////////

    //! \brief This function reorders the arrays of CFC/CFF indices and sets the pointers and sizes of the new
    //! subarrays: \details The coarse cells for interpolation from coarse to fine (iCellCFC) are divided into two
    //! subgroups: border and bulk. The fine cells (iCellCFF) are reordered accordingly. The offset cells (xOffCF,
    //! yOffCF, zOffCF) must be reordered in the same way.
    void reorderCoarseToFineIntoBorderAndBulk(int level);

    //! \brief This function reorders the arrays of FCC/FCF indices and return pointers and sizes of the new subarrays:
    //! \details The coarse cells for interpolation from fine to coarse (iCellFCC) are divided into two subgroups:
    //! border and bulk. The fine cells (iCellFCF) are reordered accordingly. The offset cells (xOffFC,
    //! yOffFC, zOffFC) must be reordered in the same way.
    void reorderFineToCoarseIntoBorderAndBulk(int level);

private:
    std::shared_ptr<GridBuilder> builder;
    std::shared_ptr<Parameter> para;
};

#endif
