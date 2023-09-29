#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../include/radial_basis.h"
#include "../../../../core/include/arrayUtils.h"
#include "../../../../core/include/vec3Operation.h"


class SwitchFuncTest : public ::testing::Test {
protected:
    double rcut;
    double rcut_smooth;
    double rji;

    static void SetUpTestSuite() {
        std::cout << "SwitchFuncTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "SwitchFuncTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 3.5;
        rcut_smooth = 3.0;
        rji = 3.5;
    }

    void TearDown() override {

    }

};  // class : SwitchFuncTest



TEST_F(SwitchFuncTest, value_rcut) {
    rcut = 3.5;
    rcut_smooth = 3.0;
    rji = 3.5;
    matersdk::mtp::SwitchFunc<double> switch_func(rcut, rcut_smooth);

    double result = switch_func.get_result(rji);
    EXPECT_DOUBLE_EQ(result, 0);
}


TEST_F(SwitchFuncTest, value_random) {
    rcut = 3.5;
    rcut_smooth = 3.0;
    rji = 3.2;
    matersdk::mtp::SwitchFunc<double> switch_func(rcut, rcut_smooth);

    double result = switch_func.get_result(rji);
    EXPECT_FLOAT_EQ(result, 0.68256);
}


TEST_F(SwitchFuncTest, value_rcut_smooth) {
    rcut = 3.5;
    rcut_smooth = 3.0;
    rji = 3.0;
    matersdk::mtp::SwitchFunc<double> switch_func(rcut, rcut_smooth);

    double result = switch_func.get_result(rji);
    EXPECT_DOUBLE_EQ(result, 1);
}




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
    matersdk::mtp::RadialBasisChebyshev<double> rb;

    rb.show_in_value();
}


TEST_F(RadialBasisChebyshevTest, constructor_1) {
    rcut = 6.0;
    rcut_smooth = 2.0;
    hmju = 3;
    rji = 3.0;
    matersdk::mtp::RadialBasisChebyshev<double> rb(rcut, rcut_smooth, hmju, rji);

    rb.show_in_value();
    rb.show_in_deriv();
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}