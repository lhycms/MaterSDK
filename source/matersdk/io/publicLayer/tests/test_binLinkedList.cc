#include <gtest/gtest.h>
#include <iostream>

#include "../include/binLinkedList.h"


class BasicStructureInfoTest : public ::testing::Test {
protected:
    int num_atoms;
    double basis_vectors[3][3];
    int atomic_numbers[12];
    double frac_coords[12][3];


    static void SetUpTestSuite() {
        std::cout << "BasicStructureInfo is setting up...\n";
    }


    static void TearDownTestSuite() {
        std::cout << "BasicStructureInfo is tearing down...\n";
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
    }

    void TearDown() override {

    }
};  // class : BasicStructureInfo 



TEST_F(BasicStructureInfoTest, default_init) {
    matersdk::BasicStructureInfo<double> basic_structure_info;
    //basic_structure_info.show();
}

TEST_F(BasicStructureInfoTest, init) {
    matersdk::Structure<double> structure_1(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::BasicStructureInfo<double> bst_1(structure_1);
    
    
    matersdk::Structure<double> structure_2;
    matersdk::BasicStructureInfo<double> bst_2(structure_2);
}


TEST_F(BasicStructureInfoTest, copy_constructor) {
    matersdk::Structure<double> structure_1(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Structure<double> structure_2;

    matersdk::BasicStructureInfo<double> bst_1(structure_1);
    matersdk::BasicStructureInfo<double> bst_2(structure_2);
    matersdk::BasicStructureInfo<double> bst_3(bst_1);
    matersdk::BasicStructureInfo<double> bst_4(bst_2);

    //bst_3.show();
    //bst_4.show();
}


TEST_F(BasicStructureInfoTest, assignment_operator) {
    matersdk::Structure<double> structure_1(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Structure<double> structure_2;

    matersdk::BasicStructureInfo<double> bst_1(structure_1);
    matersdk::BasicStructureInfo<double> bst_2(structure_2);
    matersdk::BasicStructureInfo<double> bst_3 = bst_1;
    matersdk::BasicStructureInfo<double> bst_4 = bst_2;
    bst_2 = bst_1; // or `bst_1 = bst_2;`

    //bst_3.show();
    //bst_4.show();
    //bst_1.show();
}








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



TEST_F(SupercellTest, default_constructor) {
    matersdk::Supercell<double> supercell;
    //supercell.show();
}


TEST_F(SupercellTest, constuctor_1) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    //supercell.show();
}

TEST_F(SupercellTest, assignment_operator) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    matersdk::Supercell<double> supercell_null;

    matersdk::Supercell<double> supercell_1;
    supercell_1 = supercell;
    //supercell_1.show();

    supercell_1 = supercell_null;
    //supercell_1.show();


    matersdk::Supercell<double> supercell_2(structure, scaling_matrix);
    supercell_2 = supercell;
    //supercell_2.show();
}


TEST_F(SupercellTest, copy_constructor) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    matersdk::Supercell<double> supercell;

    matersdk::Supercell<double> supercell_1(supercell);
    //supercell_1.show();
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


TEST_F(SupercellTest, get_owned_atom_idxs) {
    matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, frac_coords, false);
    scaling_matrix[0] = 3;
    scaling_matrix[1] = 3;
    scaling_matrix[2] = 1;
    matersdk::Supercell<double> supercell(structure, scaling_matrix);
    /*
        72~83   84~95   96~107
        36~47   48~59   60~71
         0~11   12~23   24~35
    */
    const int* owned_atom_idxs = supercell.get_owned_atom_idxs();
    EXPECT_EQ(owned_atom_idxs[0], 48);
    EXPECT_EQ(owned_atom_idxs[1], 49);
    EXPECT_EQ(owned_atom_idxs[2], 50);
    EXPECT_EQ(owned_atom_idxs[3], 51);
    EXPECT_EQ(owned_atom_idxs[4], 52);
    EXPECT_EQ(owned_atom_idxs[5], 53);
    EXPECT_EQ(owned_atom_idxs[6], 54);
    EXPECT_EQ(owned_atom_idxs[7], 55);
    EXPECT_EQ(owned_atom_idxs[8], 56);
    EXPECT_EQ(owned_atom_idxs[9], 57);
    EXPECT_EQ(owned_atom_idxs[10], 58);
    EXPECT_EQ(owned_atom_idxs[11], 59);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}