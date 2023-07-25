#include <gtest/gtest.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "../include/vec3Operation.h"


class Vec3OperationPointerTest : public ::testing::Test {
protected:
    double* vec1;
    double* vec2;

    static void SetUpTestSuite() {
        std::cout << "Vec3OperationTest is setting up...\n";
    }


    static void TearDownTestSuite() {
        std::cout << "Vec3OperationTest is tearing down...\n";
    }


    void SetUp() override {
        vec1 = (double*)malloc(sizeof(double) * 3);
        vec2 = (double*)malloc(sizeof(double) * 3);
        vec1[0] = 1;
        vec1[1] = 2;
        vec1[2] = 3;
        vec2[0] = 2;
        vec2[1] = 3;
        vec2[2] = 4;
    }

    
    void TearDown() override {
        free(vec1);
        free(vec2);
    }
};


TEST_F(Vec3OperationPointerTest, cross) {
    double *vertical_vec = matersdk::vec3Operation::cross(vec1, vec2);

    printf("vertical_vector = [%f, %f, %f]\n", vertical_vec[0], vertical_vec[1], vertical_vec[2]);

    free(vertical_vec);
}


TEST_F(Vec3OperationPointerTest, norm) {

}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}