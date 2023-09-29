#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../include/radial_basis.h"
#include "../../../../core/include/arrayUtils.h"
#include "../../../../core/include/vec3Operation.h"



class RadialBasisChebyshevTest : public ::testing::Test {
protected:
    double rcut;
    double rcut_smooth;
    int hmju;
    double rji;

    static void SetUpTestSuite() {
        std::cout << "RadialBasisChebyshevTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "RadialBasisChebyshevTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 6.0;
        rcut_smooth = 2.0;
        hmju = 0;
        rji = 3.0;
    }

    void TearDown() override {
    
    }
};  // class : RadialBasisChebyshevTest



TEST_F(RadialBasisChebyshevTest, constructor_default) {
    rcut = 6.0;
    rcut_smooth = 2.0;
    hmju = 3;
    rji = 3.0;
    matersdk::mtp::RadialBasisChebyshev<double> rb(rcut, rcut_smooth, hmju, rji);

    rb.show_in_value();
    rb.show_in_deriv();
}


TEST_F(RadialBasisChebyshevTest, constructor_1) {

}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}