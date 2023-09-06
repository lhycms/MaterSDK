#include "stdlib.h"

#include "lammps.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "force.h"

#include "../include/pair_se.h"


namespace matersdk {
namespace lammps {

void PairDeepPotSe::allocate() {
    // Step 1. Populate `this->allocated`
    this->allocated = 1;

    // Step 2. Populate `this->setflag`
    int num_types = this->atom->ntypes;
    this->memory->create<int>(this->setflag, num_types+1, num_types+1, "PairDeepPotSe:setflag");
    for (int ii=1; ii<=num_types; ii++) {
        for (int jj=ii; jj<=num_types; jj++)
            this->setflag[ii][jj] = 0;
    }

    // Step 3. Allocate memory
    int num_atoms = this->atom->natoms;
    int max_num_neigh_atoms = 0;
    for (int ii=0; ii<this->list->inum; ii++) {
        if (max_num_neigh_atoms < this->list->numneigh[ii])
            max_num_neigh_atoms = this->list->numneigh[ii];
    }
    this->memory->create<double>(this->relative_cart_coords, num_atoms, 3, "PairDeepPotSe:relative_cart_coord");
    this->memory->create<double>(this->tilde_r, num_atoms, max_num_neigh_atoms, 4, "PairDeepPotSe:tilde_r");
    this->memory->create<double>(this->tilde_r_deriv, num_atoms, max_num_neigh_atoms, 4, 3, "PairDeepPotSe:tilde_r_deriv");
}


// 
void PairDeepPotSe::settings(int argc, char** argv) {}


double*** PairDeepPotSe::get_neigh_relative_cart_coords() const {
    // Step 1. Allocate memory
    int num_atoms = this->atom->natoms;
    double** cart_coords = this->atom->x;
    double*** relative_cart_coords;
    int num_max_neigh_atoms = 0;
    for (int ii=0; ii<this->list->inum; ii++) {
        if (num_max_neigh_atoms < this->list->numneigh[ii]) 
            num_max_neigh_atoms = this->list->numneigh[ii];
    }
    this->memory->create<double>(relative_cart_coords, this->list->inum, num_max_neigh_atoms, 3);

    int center_index, neigh_index;
    double* center_cart_coord;
    this->memory->create<double>(center_cart_coord, 3, "PairDeepPotSe:center_cart_coord");
    double* neigh_cart_coord;
    this->memory->create<double>(neigh_cart_coord, 3, "PairDeepPotSe:neigh_cart_coord");

    // Step 2. Populate `relative_cart_coords`
    for (int ii=0; ii<this->list->inum; ii++) {
        center_index = this->list->ilist[ii];
        center_cart_coord[0] = cart_coords[center_index][0];
        center_cart_coord[1] = cart_coords[center_index][1];
        center_cart_coord[2] = cart_coords[center_index][2];
        for (int jj=0; jj<this->list->numneigh[ii]; jj++) {
            neigh_index = this->list->firstneigh[ii][jj];
            neigh_cart_coord[0] = cart_coords[neigh_index][0];
            neigh_cart_coord[1] = cart_coords[neigh_index][1];
            neigh_cart_coord[2] = cart_coords[neigh_index][2];
            relative_cart_coords[ii][jj][0] = neigh_cart_coord[0] - center_cart_coord[0];
            relative_cart_coords[ii][jj][1] = neigh_cart_coord[1] - center_cart_coord[1];
            relative_cart_coords[ii][jj][2] = neigh_cart_coord[2] - center_cart_coord[2];
        }
    }

    // Step . Free memory
    this->memory->destroy<double>(center_cart_coord);
    this->memory->destroy<double>(neigh_cart_coord);

    return relative_cart_coords;
}



}   // namespace: lammps
}   // namespace: matersdk 