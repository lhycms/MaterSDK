#include <gtest/gtest.h>
#include <iostream>
#include <stdlib.h>

#include "../include/structure.h"


class StructureTest : public ::testing::Test {
protected:
    matersdk::Structure<double> *ptr_structure;


    static void SetUpTestSuite() {
        std::cout << "StructureTest is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "StructureTest is tearing down...\n";
    }

    void SetUp() override {
        int num_atoms = 12;
        double basis_vectors[3][3] = {
            {3.1903157348, 5.5257885468, 0.0000000000},
            {-6.3806307800, 0.0000000000, 0.0000000000},
            {0.0000000000, 0.0000000000, 23.1297687334}
        };
        int atomic_numbers[12] = {
                    42, 16, 16, 42,
                    16, 16, 42, 16,
                    16, 42, 16, 16
        };
        double cart_coords[12][3] = {
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0} ,
            {1.0, 2.0, 3.0}
        };


        matersdk::Structure<double> structure(num_atoms, basis_vectors, atomic_numbers, cart_coords);
    }

    void TearDown() override {
        
    }
};



TEST_F(StructureTest, init) {
    std::cout << 111 << std::endl;
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}