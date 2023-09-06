#include "stdlib.h"

#include "lammps.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "force.h"

#include "torch/torch.h"
#include "torch/script.h"

#include "./pair_deeppot_se.h"



namespace LAMMPS_NS {


PairDeepPotSe::PairDeepPotSe(LAMMPS* lmp) : Pair(lmp) {

}


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
    
    this->relative_cart_coords = this->get_neigh_relative_cart_coords();
}


// Parsing the parameter after `pair_style` which is `global potential parameters`
/*
    e.g. pair_style deeppot/se 6.5
            - deeppot/se: 势函数类型
            - 6.5: 截断半径
*/
void PairDeepPotSe::settings(int argc, char** argv) {
    if (argc != 1)
        this->error->all(FLERR, "Illegal pair_style command");
    
    this->rcut = utils::numeric(FLERR, argv[0], false, lmp);
}


// Parsing the `coefficients of pair potential`
/*
    e.g. pair_coeff * * model.pt
*/
void PairDeepPotSe::coeff(int argc, char** argv) {
    int ilo, ihi, hlo, jhi;
    utils::bounds(FLERR, argv[0], 1, this->atom->ntypes, ilo, ihi, this->error);
    utils::bounds(FLERR, argv[1], 1, this->atom->ntypes, jlo, jhi, this->error);

    // Step 1. Load the NN model
    try {
        this->model = torch::jit::load(argv[2]);
    } catch (const c10::Error& e) {
        this->error->all(FLERR, "Error loading the DeepPot-SE model");
    }

    // Step 2. Read the rcut for all atom types. 
    //      Note. `this->setflag[ii][ii]` must be set. -> `this->setflag[ii][ii] = 1;`
    int num_type_pairs;
    for (int ii=ilo; ii<=ihi; ii++) {
        for (int jj= MAX(jlo, ii); jj<=jhi; jj++) {
            this->setflag[ii][jj] = 1;
            this->cut[ii][jj] = this->rcut;
            this->cutsq[ii][jj] = this->rcut * this->rcut;
            num_type_pairs++;
        }
    }

    if (num_type_pairs == 0)
        this->error->all(FLERR, "Incorrect args for pair coefficients");
}


// Assign pairwise coefficients for all atom types
// Init for one type pair i,j and corresponding j, i.
double PairDeepPotSe::init_one(int ii, int jj) {
    if (this->setflag[ii][jj] == 0)
        this->error->all(FLERR, "All pair coeffs are not set");
    
    this->cut[jj][ii] = this->cut[ii][jj];

    return this->cut[ii][jj];
}


// add_request a NeighList
void PairDeepPotSe::init_style() {
    if (this->force->newton_pair == 0)
        this->error->all(FLERR, "Pair style deeppot/se requires newton_pair on");

    this->neighbor->add_request(this, NeighConst::REQ_FULL);
}


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


}   // namespace: LAMMPS_NS