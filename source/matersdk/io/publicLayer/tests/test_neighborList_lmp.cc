#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>

#include "../include/structure.h"
#include "../include/neighborList.h"
#include "../../../../core/include/vec3Operation.h"
#include "../../../../core/include/arrayUtils.h"


class NeighborListTest : public ::testing::Test {
protected:
    int num_atoms;
    double basis_vectors[3][3];
    int atomic_numbers[12];
    double frac_coords[12][3];
    double rcut;
    //double bin_size_xyz[3];
    bool pbc_xyz[3];

    matersdk::Structure<double> structure;
    matersdk::NeighborList<double> neighbor_list;

    // Supercell properties
    int inum;
    int* ilist;
    int* numneigh;
    int** firstneigh;
    int* types;
    double** x;

    // 
    int num_center_atomic_numbers;
    int* center_atomic_numbers_lst;
    int num_neigh_atomic_numbers;
    int* neigh_atomic_numbers_lst;
    int* num_neigh_atoms_lst;
    int tot_num_neigh_atoms;


    static void SetUpTestSuite() {
        std::cout << "NeighborListTest is setting up...\n";
    }


    static void TearDownTestSuite() {
        std::cout << "NeighborListTest is tearing down...\n";
    }


    void SetUp() override {
        num_atoms = 12;

        basis_vectors[0][0] = 3.1903157348;
        basis_vectors[0][1] = 5.5257885468;
        basis_vectors[0][2] = 0.0000000000;
        basis_vectors[1][0] = -6.3806307800;
        basis_vectors[1][1] = 0.0000000000;
        basis_vectors[1][2] = 0.0000000000;
        basis_vectors[2][0] = 0.0000000000;
        basis_vectors[2][1] = 0.0000000000;
        basis_vectors[2][2] = 23.1297687334;

        atomic_numbers[0] = 42;
        atomic_numbers[1] = 16;
        atomic_numbers[2] = 16;
        atomic_numbers[3] = 42;
        atomic_numbers[4] = 16;
        atomic_numbers[5] = 16;
        atomic_numbers[6] = 42;
        atomic_numbers[7] = 16;
        atomic_numbers[8] = 16;
        atomic_numbers[9] = 42; 
        atomic_numbers[10] = 16;
        atomic_numbers[11] = 16;

        frac_coords[0][0] = 0.333333333333;
        frac_coords[0][1] = 0.166666666667;
        frac_coords[0][2] = 0.500000000000;
        frac_coords[1][0] = 0.166666666667;
        frac_coords[1][1] = 0.333333333333;
        frac_coords[1][2] = 0.432343276548;
        frac_coords[2][0] = 0.166666666667;
        frac_coords[2][1] = 0.333333333333;
        frac_coords[2][2] = 0.567656723452;
        frac_coords[3][0] = 0.333333333333;
        frac_coords[3][1] = 0.666666666667;
        frac_coords[3][2] = 0.500000000000;
        frac_coords[4][0] = 0.166666666667;
        frac_coords[4][1] = 0.833333333333;
        frac_coords[4][2] = 0.432343276548;
        frac_coords[5][0] = 0.166666666667;
        frac_coords[5][1] = 0.833333333333;
        frac_coords[5][2] = 0.567656723452;
        frac_coords[6][0] = 0.833333333333;
        frac_coords[6][1] = 0.166666666667;
        frac_coords[6][2] = 0.500000000000;
        frac_coords[7][0] = 0.666666666667;
        frac_coords[7][1] = 0.333333333333;
        frac_coords[7][2] = 0.432343276548;
        frac_coords[8][0] = 0.666666666667;
        frac_coords[8][1] = 0.333333333333;
        frac_coords[8][2] = 0.567656723452;
        frac_coords[9][0] = 0.833333333333;
        frac_coords[9][1] = 0.666666666667;
        frac_coords[9][2] = 0.500000000000;
        frac_coords[10][0] = 0.666666666667;
        frac_coords[10][1] = 0.833333333333;
        frac_coords[10][2] = 0.432343276548;
        frac_coords[11][0] = 0.666666666667;
        frac_coords[11][1] = 0.833333333333;
        frac_coords[11][2] = 0.567656723452;

        rcut = 3.3;
        pbc_xyz[0] = true;
        pbc_xyz[1] = true;
        pbc_xyz[2] = false; 


        structure = matersdk::Structure<double>(
                        num_atoms, 
                        basis_vectors, 
                        atomic_numbers, 
                        frac_coords, 
                        false);
        neighbor_list = matersdk::NeighborList<double>(
                        structure,
                        rcut,
                        pbc_xyz,
                        true);

        // Supercell properties simulating lammps
        inum = neighbor_list.get_binLinkedList().get_supercell().get_prim_num_atoms();
        
        ilist = (int*)malloc(sizeof(int) * inum);
        for (int ii=0; ii<inum; ii++)
            ilist[ii] = 0;
        for (int ii=0; ii<inum; ii++) 
            ilist[ii] = ii;

        numneigh = (int*)malloc(sizeof(int) * inum);
        for (int ii=0; ii<inum; ii++) 
            numneigh[ii] = 0;
        for (int ii=0; ii<inum; ii++)
            numneigh[ii] = neighbor_list.get_neighbor_lists()[ii].size();

        firstneigh = (int**)malloc(sizeof(int*) * inum);
        for (int ii=0; ii<inum; ii++)
            firstneigh[ii] = (int*)malloc(sizeof(int) * numneigh[ii]);
        for (int ii=0; ii<inum; ii++) {
            for (int jj=0; jj<numneigh[ii]; jj++) {
                firstneigh[ii][jj] = 0;
            }
        }
        for (int ii=0; ii<inum; ii++) {
            for (int jj=0; jj<numneigh[ii]; jj++)
                firstneigh[ii][jj] = neighbor_list.get_neighbor_lists()[ii][jj];
        }


        types = (int*)neighbor_list.get_binLinkedList().get_supercell().get_structure().get_atomic_numbers();
        x = (double**)neighbor_list.get_binLinkedList().get_supercell().get_structure().get_cart_coords();

        //
        num_center_atomic_numbers = 2;
        center_atomic_numbers_lst = (int*)malloc(sizeof(int) * num_neigh_atomic_numbers);
        center_atomic_numbers_lst[0] = 16;
        center_atomic_numbers_lst[1] = 42;
        num_neigh_atomic_numbers = 2;
        neigh_atomic_numbers_lst = (int*)malloc(sizeof(int) * num_neigh_atomic_numbers);
        neigh_atomic_numbers_lst[0] = 16;
        neigh_atomic_numbers_lst[1] = 42;
        num_neigh_atoms_lst = (int*)malloc(sizeof(int) * num_neigh_atomic_numbers);
        num_neigh_atoms_lst[0] = 8;
        num_neigh_atoms_lst[1] = 7;
        tot_num_neigh_atoms = 0;
        for (int ii=0; ii<num_neigh_atomic_numbers; ii++) {
            tot_num_neigh_atoms += num_neigh_atoms_lst[ii];
        }
    }


    void TearDown() override {
        free(ilist);
        free(numneigh);
        for (int ii=0; ii<inum; ii++)
            free(firstneigh[ii]);
        free(firstneigh);
        free(center_atomic_numbers_lst);
        free(neigh_atomic_numbers_lst);
        free(num_neigh_atoms_lst);
    }

}; // class: NeighborListTest


TEST_F(NeighborListTest, assign_list_neigh4fortran) {
    // Step 1. Allocate memory for `list_neigh`
    int** list_neigh = (int**)malloc(sizeof(int*) * inum);
    for (int ii=0; ii<inum; ii++) {
        list_neigh[ii] = (int*)malloc(sizeof(int) * tot_num_neigh_atoms);
    }
    for (int ii=0; ii<inum; ii++) {
        for (int jj=0; jj<tot_num_neigh_atoms; jj++)
            list_neigh[ii][jj] = 0;
    }

    // Step 2. Populate `list_neigh`
    matersdk::NeighborList<double>::assign_list_neigh4fortran(
                list_neigh,
                inum,
                ilist,
                numneigh,
                firstneigh,
                types,
                x,
                num_center_atomic_numbers,
                center_atomic_numbers_lst,
                num_neigh_atomic_numbers,
                neigh_atomic_numbers_lst,
                num_neigh_atoms_lst);

    // Step 3. Print out `type` catched by `list_neigh`
    int* sorted_center_an = (int*)malloc(sizeof(int) * inum);   // 16, 16, 16, 16, 16, 16, 16, 16, 42, 42, 42, 42
    int tmp_center_an;
    int tmp_cidx = 0;
    for (int ii=0; ii<num_center_atomic_numbers; ii++) {
        for (int jj=0; jj<inum; jj++) {
            tmp_center_an = ilist[jj];
            if (types[tmp_center_an] != center_atomic_numbers_lst[ii])
                continue;
            sorted_center_an[tmp_cidx] = types[tmp_center_an];
            tmp_cidx++;
        }
    }

    // Step 3.1. 
    for (int ii=0; ii<inum; ii++) {
        printf("%4d : ", sorted_center_an[ii]);
        for (int jj=0; jj<tot_num_neigh_atoms; jj++) {
            if (list_neigh[ii][jj] != 0)
                printf("%3d, ", types[list_neigh[ii][jj]-1]);
            else
                printf("%3d, ", 0);
        }
        printf("\n");
    }
    
    // Step . Free memory
    for (int ii=0; ii<inum; ii++)
        free(list_neigh[ii]);
    free(list_neigh);
    free(sorted_center_an);
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}