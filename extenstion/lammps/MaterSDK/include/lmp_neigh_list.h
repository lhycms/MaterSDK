#ifndef MATERSDK_LMP_SE_H
#define MATERSDK_LMP_SE_H

#include "lammps.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "memory.h"


namespace matersdk {
namespace lammps {

class LmpNeighList : public LAMMPS_NS::Pair{
public:
    // Return.shape = (num_atoms, num_nbrs, 3)
    double*** get_neigh_relative_cart_coords() const;

private:
    LAMMPS_NS::LAMMPS* lmp;
};  // class: LmpNeighList


}   // namespace: lammps
}   // namespace: matersdk

#endif