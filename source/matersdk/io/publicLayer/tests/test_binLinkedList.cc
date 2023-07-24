#include <gtest/gtest.h>
#include <iostream>

#include "../include/binLinkedList.h"



class SupercellTest : public ::testing::Test {
protected:
    int num_atoms;
    double basis_vectors[3][3];
    int atomic_numbers[12];
    double frac_coords[12][3];
    int scaling_matrix[3];

    static void SetUpTestSuite() {
        std::cout << "SupercellTest is setting up...\n";
    }


    static void TearDownTestSuite() {
        std::cout << "SupercellTest is tearing down...\n";
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

        scaling_matrix[0] = 3;
        scaling_matrix[1] = 3;
        scaling_matrix[2] = 1;
    }


    void TearDown() override {
        
    }
}; // class: Supercell class



TEST_F(SupercellTest, init) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
}


TEST_F(SupercellTest, calc_prim_cell_idx_xyz_even) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    // Test 1: Scaling factor 是奇数
    scaling_matrix[0] = 5;
    scaling_matrix[1] = 7;
    scaling_matrix[2] = 9;
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    // supercell.calc_prim_cell_idx_xyz();
    // supercell.calc_prim_cell_idx();
    const int* prim_cell_idx_xyz = supercell.get_prim_cell_idx_xyz();

    EXPECT_EQ(prim_cell_idx_xyz[0], 2);
    EXPECT_EQ(prim_cell_idx_xyz[1], 3);
    EXPECT_EQ(prim_cell_idx_xyz[2], 4);

    const int prim_cell_idx = supercell.get_prim_cell_idx();
    std::cout << prim_cell_idx << std::endl;
}


TEST_F(SupercellTest, calc_cell_idx_xyz_odd) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    // Test 1: Scaling factor 是偶数
    scaling_matrix[0] = 6;
    scaling_matrix[1] = 8;
    scaling_matrix[2] = 10;
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    // supercell.calc_prim_cell_idx_xyz();
    // supercell.calc_prim_cell_idx();
    const int* prim_cell_idx_xyz = supercell.get_prim_cell_idx_xyz();
    EXPECT_EQ(prim_cell_idx_xyz[0], 2);
    EXPECT_EQ(prim_cell_idx_xyz[1], 3);
    EXPECT_EQ(prim_cell_idx_xyz[2], 4);

    const int prim_cell_idx = supercell.get_prim_cell_idx();
    std::cout << prim_cell_idx << std::endl;
}


TEST_F(SupercellTest, get_num_atoms) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    scaling_matrix[0] = 6;
    scaling_matrix[1] = 8;
    scaling_matrix[2] = 10;
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    
    EXPECT_EQ(supercell.get_prim_num_atoms(), 12);
    EXPECT_EQ(supercell.get_num_atoms(), 12 * scaling_matrix[0] * scaling_matrix[1] * scaling_matrix[2]);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}