#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>

#include "../../../io/publicLayer/include/structure.h"
#include "../include/se.h"


class PairTildeRTest : public ::testing::Test {
protected:
    int num_atoms;
    double basis_vectors[3][3];
    int atomic_numbers[12];
    double frac_coords[12][3];
    double rcut;
    double bin_size_xyz[3];
    bool pbc_xyz[3];  

    int center_atomic_number;
    int neigh_atomic_number;
    int num_neigh_atoms;
    double rcut_smooth;

    matersdk::Structure<double> structure;
    matersdk::NeighborList<double> neighbor_list;


    static void SetUpTestSuite() {
        std::cout << "PairTildeR TestSuite is setting up...\n";
    }


    static void TearDownTestSuite() {
        std::cout << "PairTildeR TestSuite is tearing down...\n";
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
        bin_size_xyz[0] = 3.0;
        bin_size_xyz[1] = 3.0;
        bin_size_xyz[2] = 3.0;
        pbc_xyz[0] = true;
        pbc_xyz[1] = true;
        pbc_xyz[2] = false; 


        structure = matersdk::Structure<double>(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
        neighbor_list = matersdk::NeighborList<double>(structure, rcut, pbc_xyz, true);
    }

    void TearDown() override {

    }
};


TEST_F(PairTildeRTest, default_constructor) {
    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r;
    //pair_tilde_r.show();
}


TEST_F(PairTildeRTest, constructor_1) {
    center_atomic_number = 42;
    neigh_atomic_number = 16;
    num_neigh_atoms = 24;
    rcut_smooth = 3.0;
    
    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r(neighbor_list, center_atomic_number, neigh_atomic_number, num_neigh_atoms, rcut_smooth);
    //pair_tilde_r.show();
}


TEST_F(PairTildeRTest, constructor_2) {
    center_atomic_number = 42;
    neigh_atomic_number = 16;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r(neighbor_list, center_atomic_number, neigh_atomic_number, rcut_smooth);
    //pair_tilde_r.show();
}


TEST_F(PairTildeRTest, constructor_3) {
    pbc_xyz[0] = true;
    pbc_xyz[1] = true;
    pbc_xyz[2] = false;
    bool sort = true;

    center_atomic_number = 42;
    neigh_atomic_number = 16;
    num_neigh_atoms = 100;
    rcut = 3.3;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r(structure, rcut, pbc_xyz, sort, center_atomic_number, neigh_atomic_number, num_neigh_atoms, rcut_smooth);
    //pair_tilde_r.show();
}


TEST_F(PairTildeRTest, constructor_4) {
    pbc_xyz[0] = true;
    pbc_xyz[1] = true;
    pbc_xyz[2] = false;
    bool sort = true;

    center_atomic_number = 16;
    neigh_atomic_number = 16;
    rcut = 3.3;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r(structure, rcut, pbc_xyz, sort, center_atomic_number, neigh_atomic_number, rcut_smooth);
    pair_tilde_r.show();
}


TEST_F(PairTildeRTest, get_num_center_atoms) {
    // Case 1.
    center_atomic_number = 42;
    neigh_atomic_number = 16;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r_42_16(neighbor_list, center_atomic_number, neigh_atomic_number, rcut_smooth);
    double num_center_atoms_42_16 = pair_tilde_r_42_16.get_num_center_atoms();
    EXPECT_EQ(num_center_atoms_42_16, 4);

    // Case 2.
    center_atomic_number = 16;
    neigh_atomic_number = 16;
    rcut_smooth = 3.0;
    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r_16_16(neighbor_list, center_atomic_number, neigh_atomic_number, rcut_smooth);
    double num_center_atoms_16_16 = pair_tilde_r_16_16.get_num_center_atoms();
    EXPECT_EQ(num_center_atoms_16_16, 8);
}


TEST_F(PairTildeRTest, get_num_neigh_atoms) {
    // Case 1.
    center_atomic_number = 16;
    neigh_atomic_number = 16;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r_16_16(neighbor_list, center_atomic_number, neigh_atomic_number, rcut_smooth);
    double num_neigh_atoms_16_16 = pair_tilde_r_16_16.get_num_neigh_atoms();
    EXPECT_EQ(num_neigh_atoms_16_16, 7);

    // Case 2.
    center_atomic_number = 16;
    neigh_atomic_number = 16;
    num_neigh_atoms = 100;
    rcut_smooth = 3.0;

    matersdk::deepPotSE::PairTildeR<double> pair_tilde_r_16_16_(neighbor_list, center_atomic_number, neigh_atomic_number, num_neigh_atoms, rcut_smooth);
    double num_neigh_atoms_16_16_ = pair_tilde_r_16_16_.get_num_neigh_atoms();
    EXPECT_EQ(num_neigh_atoms_16_16_, 100);

}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}