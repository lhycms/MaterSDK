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

class LmpDeepPotSe : public LAMMPS_NS::Pair{
public:
    double** relative_cart_coords;
    double*** tilde_r;
    double**** tilde_r_deriv;

    void allocate();

    void settings(int argc, char** argv);

    void deallocate();

    // Return.shape = (num_atoms, num_nbrs, 3)
    double*** get_neigh_relative_cart_coords() const;

    // 
    
};  // class: LmpDeepPotSe


}   // namespace: lammps
}   // namespace: matersdk

#endif