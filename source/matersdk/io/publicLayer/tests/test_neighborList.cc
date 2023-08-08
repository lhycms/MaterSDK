#include <gtest/gtest.h>
#include <iostream>

#include "../include/structure.h"
#include "../include/neighborList.h"
#include "../../../../core/include/vec3Operation.h"


class NeighborListTest : public ::testing::Test {
protected:
    int num_atoms;
    double basis_vectors[3][3];
    int atomic_numbers[12];
    double frac_coords[12][3];
    double rcut;
    double bin_size_xyz[3];
    bool pbc_xyz[3];


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
        bin_size_xyz[0] = 3.0;
        bin_size_xyz[1] = 3.0;
        bin_size_xyz[2] = 3.0;
        pbc_xyz[0] = true;
        pbc_xyz[1] = true;
        pbc_xyz[2] = false; 
    }


    void TearDown() override {

    }
};


TEST_F(NeighborListTest, constructor_1) {
    rcut = 3.3;             // 截断半径
    bin_size_xyz[0] = 1.65;  // X 方向上的 bin_size (一般默认 rcut/2)
    bin_size_xyz[1] = 1.65;  // Y 方向上的 bin_size
    bin_size_xyz[2] = 1.65;  // Z 方向上的 bin_size
    pbc_xyz[0] = true;      // X 方向上是否具有周期性
    pbc_xyz[1] = true;      // Y 方向上是否具有周期性
    pbc_xyz[2] = false;     // Z 方向上是否具有周期性
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, bin_size_xyz, pbc_xyz, true);
    
    //neighbor_list.show_in_index();
    //printf("\n");
    //neighbor_list.show_in_prim_index();
    //printf("\n");
    neighbor_list.show_in_an();
    printf("\n");
    neighbor_list.show_in_distances();
}


TEST_F(NeighborListTest, constructor_2) {
    rcut = 3.3;             // 截断半径
    pbc_xyz[0] = true;      // X 方向上是否具有周期性
    pbc_xyz[1] = true;      // Y 方向上是否具有周期性
    pbc_xyz[2] = false;     // Z 方向上是否具有周期性
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, pbc_xyz, true);
    
    //neighbor_list.show_in_index();
    //printf("\n");
    //neighbor_list.show_in_prim_index();
    //printf("\n");
    neighbor_list.show_in_an();
    printf("\n");
    neighbor_list.show_in_distances();
}


TEST_F(NeighborListTest, get_max_num_neigh_atoms) {
    rcut = 3.3;             // 截断半径
    bin_size_xyz[0] = 3.0;  // X 方向上的 bin_size (一般默认 rcut/2)
    bin_size_xyz[1] = 3.0;  // Y 方向上的 bin_size
    bin_size_xyz[2] = 3.0;  // Z 方向上的 bin_size
    pbc_xyz[0] = true;      // X 方向上是否具有周期性
    pbc_xyz[1] = true;      // Y 方向上是否具有周期性
    pbc_xyz[2] = false;     // Z 方向上是否具有周期性
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, bin_size_xyz, pbc_xyz);

    const int max_num_neigh_atoms = neighbor_list.get_max_num_neigh_atoms();
    EXPECT_EQ(max_num_neigh_atoms, 12);
}


TEST_F(NeighborListTest, get_max_num_neigh_atoms_ss) {
    rcut = 2.8;
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, bin_size_xyz, pbc_xyz);

    const int max_num_neigh_atoms_42 = neighbor_list.get_max_num_neigh_atoms_ss(42);
    const int max_num_neigh_atoms_16 = neighbor_list.get_max_num_neigh_atoms_ss(16);

    EXPECT_EQ(max_num_neigh_atoms_42, 3);
    EXPECT_EQ(max_num_neigh_atoms_16, 6);
}


TEST_F(NeighborListTest, get_num_neigh_atoms) {
    rcut = 3.3;
    pbc_xyz[0] = true;
    pbc_xyz[1] = true;
    pbc_xyz[2] = false;

    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, pbc_xyz);
    const int num_neigh_atoms_42 = neighbor_list.get_num_neigh_atoms(0);
    const int num_neigh_atoms_16 = neighbor_list.get_num_neigh_atoms(1);
    
    EXPECT_EQ(num_neigh_atoms_42, 12);
    EXPECT_EQ(num_neigh_atoms_16, 10);
}


TEST_F(NeighborListTest, get_neigh_atomic_numbers) {
    rcut = 3.3;
    pbc_xyz[0] = true;      
    pbc_xyz[1] = true;      
    pbc_xyz[2] = false; 
    
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, pbc_xyz, true);

    // Step 1. Run the program
    int* neigh_atomic_numbers_42 = neighbor_list.get_neigh_atomic_numbers(0);
    int* neigh_atomic_numbers_16 = neighbor_list.get_neigh_atomic_numbers(1);

    // Step 2. calculate the standard value
    int standard_neigh_atomic_numbers_42[12] = { 16, 16, 16, 16, 16, 16, 42, 42, 42, 42, 42, 42 };
    int standard_neigh_atomic_numbers_16[10] = { 42, 42, 42, 16, 16, 16, 16, 16, 16, 16 };

    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(0); ii++) {
        EXPECT_EQ(neigh_atomic_numbers_42[ii], standard_neigh_atomic_numbers_42[ii]);
    }
    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(1); ii++) {
        EXPECT_EQ(neigh_atomic_numbers_16[ii], standard_neigh_atomic_numbers_16[ii]);
    }

    // Step . Free memory
    free(neigh_atomic_numbers_42);
    free(neigh_atomic_numbers_16);
}


TEST_F(NeighborListTest, get_neigh_distances) {
    rcut = 3.3;
    pbc_xyz[0] = true;
    pbc_xyz[1] = true;
    pbc_xyz[2] = false;

    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, pbc_xyz, true);

    // Step 1. Run the program
    double* distances_42 = neighbor_list.get_neigh_distances(0);
    double* distances_16 = neighbor_list.get_neigh_distances(1);

    // Step 2. Calculate the standard value
    double standard_distances_42[] = 
        { 2.416933, 2.416933, 2.416933, 2.416933, 2.416933, 2.416933, 3.190315, 3.190315, 3.190315, 3.190315, 3.190316, 3.190316 };
    double standard_distances_16[] = 
        { 2.416933, 2.416933, 2.416933, 3.129769, 3.190315, 3.190315, 3.190315, 3.190315, 3.190316, 3.190316 };

    // Step Test
    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(0); ii++) {
        EXPECT_FLOAT_EQ(distances_42[ii], standard_distances_42[ii]);
    }
    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(1); ii++) {
        EXPECT_FLOAT_EQ(distances_16[ii], standard_distances_16[ii]);
    }

    // Step . Free memory
    free(distances_42);
    free(distances_16);
}


TEST_F(NeighborListTest, get_neigh_relative_cart_coords) {
    rcut = 3.3;
    pbc_xyz[0] = true;
    pbc_xyz[1] = true;
    pbc_xyz[2] = false;

    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, pbc_xyz, true);

    double** relative_cart_coords = neighbor_list.get_neigh_relative_cart_coords(0);
    double tmp_distance;

    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(0); ii++) {
        tmp_distance = matersdk::vec3Operation::norm(relative_cart_coords[ii]);
        printf("%f, ", tmp_distance);
    }
    printf("\n");


    // Step . Free memory
    for (int ii=0; ii<neighbor_list.get_num_neigh_atoms(0); ii++) {
        free(relative_cart_coords[ii]);
    }
    free(relative_cart_coords);

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}