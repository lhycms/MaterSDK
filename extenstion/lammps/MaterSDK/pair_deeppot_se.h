#ifndef MATERSDK_LMP_SE_H
#define MATERSDK_LMP_SE_H

#include "lammps.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "memory.h"

#include <torch/torch.h>
#include <torch/script.h>


namespace LAMMPS_NS {


class PairDeepPotSe : public LAMMPS_NS::Pair{
public:
    double rcut;    // 截断半径
    double** cut;   // 截断半径
    double** relative_cart_coords;
    double*** tilde_r;
    double**** tilde_r_deriv;
    torch::jit::script::Module model;

    // Constructor
    PairDeepPotSe(LAMMPS_NS::LAMMPS* lmp);

    // Allocates memory to arrays used to calculate pair interaction forces.
    void allocate();

    // Reads and process `global pair potential parameters` entered after `pair_style` command in LAMMPS input script.
    void settings(int argc, char** argv);

    // Reads and process the coefficients of `pair potential parameter`
    void coeff(int argc, char** argv);

    // Assigns pairwise coefficients for all atom types
    double init_one(int i, int j);

    // Request NeighList
    void init_style();

    // 
    void deallocate();

    // Return.shape = (num_atoms, num_nbrs, 3)
    double*** get_neigh_relative_cart_coords() const;

    // 
    
};  // class: LmpDeepPotSe



}   // namespace: LAMMPS_NS

#endif