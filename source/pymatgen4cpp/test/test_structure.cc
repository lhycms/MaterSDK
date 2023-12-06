#include <gtest/gtest.h>
#include <iostream>
#include <stdlib.h>

#include "../include/load_file.h"



class StructureTest : public ::testing::Test {
protected:    
    static void SetUpTestSuite() {
        std::cout << "StructureTest is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "StructureTest is tearing down...\n";
    }

    void SetUp() override {
    }

    void TearDown() override {
    }
}; // class: StructureTest



TEST_F(StructureTest, from_file) {    
    PymatgenUtils::Structure structure;
    structure.from_file("/data/home/liuhanyu/hyliu/pwmat_demo/MoS2/scf_/POSCAR");
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}