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


TEST_F(NeighborListTest, copy_constructor) {
    rcut = 3.3;           
    pbc_xyz[0] = true;    
    pbc_xyz[1] = true;     
    pbc_xyz[2] = false;
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    
    matersdk::NeighborList<double> neighbor_list_1;
    matersdk::NeighborList<double> neighbor_list_2(structure, rcut, pbc_xyz, false);

    matersdk::NeighborList<double> neighbor_list_3(neighbor_list_1);
    matersdk::NeighborList<double> neighbor_list_4(neighbor_list_2);

    neighbor_list_3.show_in_an();
    neighbor_list_4.show_in_an();
}


TEST_F(NeighborListTest, assignment_operator) {
    rcut = 3.3;           
    pbc_xyz[0] = true;    
    pbc_xyz[1] = true;     
    pbc_xyz[2] = false;
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);

    matersdk::NeighborList<double> neighbor_list_1;
    matersdk::NeighborList<double> neighbor_list_2(structure, rcut, pbc_xyz, false);
    matersdk::NeighborList<double> neighbor_list_3;
    matersdk::NeighborList<double> neighbor_list_4(structure, rcut, pbc_xyz, false);

    neighbor_list_3 = neighbor_list_1;
    //neighbor_list_3 = neighbor_list_2;
    neighbor_list_4 = neighbor_list_1;
    //neighbor_list_4 = neighbor_list_2;
    
    neighbor_list_4.show_in_an();
}


TEST_F(NeighborListTest, get_num_center_atoms) {
    rcut = 3.3;             // 截断半径
    bin_size_xyz[0] = 3.0;  // X 方向上的 bin_size (一般默认 rcut/2)
    bin_size_xyz[1] = 3.0;  // Y 方向上的 bin_size
    bin_size_xyz[2] = 3.0;  // Z 方向上的 bin_size
    pbc_xyz[0] = true;      // X 方向上是否具有周期性
    pbc_xyz[1] = true;      // Y 方向上是否具有周期性
    pbc_xyz[2] = false;     // Z 方向上是否具有周期性
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, bin_size_xyz, pbc_xyz);

    const int num_center_atoms = neighbor_list.get_num_center_atoms();
    EXPECT_EQ(num_center_atoms, 12);
}


TEST_F(NeighborListTest, get_rcut) {
    rcut = 3.3;             // 截断半径
    bin_size_xyz[0] = 3.0;  // X 方向上的 bin_size (一般默认 rcut/2)
    bin_size_xyz[1] = 3.0;  // Y 方向上的 bin_size
    bin_size_xyz[2] = 3.0;  // Z 方向上的 bin_size
    pbc_xyz[0] = true;      // X 方向上是否具有周期性
    pbc_xyz[1] = true;      // Y 方向上是否具有周期性
    pbc_xyz[2] = false;     // Z 方向上是否具有周期性
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::NeighborList<double> neighbor_list(structure, rcut, bin_size_xyz, pbc_xyz);

    const double rcut = neighbor_list.get_rcut();
    EXPECT_FLOAT_EQ(rcut, 3.3);
}






int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}