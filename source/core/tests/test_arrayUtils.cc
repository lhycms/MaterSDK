#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>

#include "../include/arrayUtils.h"



class ArrayUtilsTest : public ::testing::Test {
protected:
    int element_size_0;
    int element_size_1;
    int element_size_2;
    bool init_mark;

    static void SetUpTestSuite() {
        std::cout << "ArrayUtilsTest TestSuite is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "ArrayUtilsTest TestSuite is tearing down...\n";
    }

    void SetUp() override  {
        element_size_0 = 3;
        element_size_1 = 4;
        element_size_2 = 5;
        init_mark = true;
    }

    void TearDown() override {

    }
}; // class : ArrayUitlsTest



TEST_F(ArrayUtilsTest, array_3d) {
    init_mark = true;
    
    double*** pointer_3dArray = matersdk::arrayUtils::allocate3dArray<double>(element_size_0, element_size_1, element_size_2, init_mark);
    matersdk::arrayUtils::show3dArray(pointer_3dArray, element_size_0, element_size_1, element_size_2);
    matersdk::arrayUtils::free3dArray(pointer_3dArray, element_size_0, element_size_1);
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}