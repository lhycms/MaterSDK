#include <gtest/gtest.h>
#include <iostream>

// cmake -DBUILD_TEST=1 ..; make -j 8; ./bin/core/test_AlignedArray
#include "../include/AlignedArray.h"


class AlignedArrayTest : public ::testing::Test {
protected:
    matersdk::AlignedArray<double> *ptr_aa;

    static void SetUpTestSuite() {
        std::cout << "Set up the AlignedArrayTest (TestSuite)...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "Tear down the AlignedArragTest (TestSuite)...\n";
    }

    void SetUp() override {
        ptr_aa = new matersdk::AlignedArray<double>(100);
    }


    void TearDown() override {
        delete ptr_aa;
    }
};


TEST_F(AlignedArrayTest, Test_1) {
    std::cout << ptr_aa->size() << std::endl;
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}