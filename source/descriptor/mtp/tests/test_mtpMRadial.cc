#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include "../include/mtpMRadial.h"


class MtpSwitchFunc1Test : public ::testing::Test
{
protected:
    double rcut;
    double rcut_smooth;
    double distance_ij;    // absolute value
    static void SetUpTestSuite() {
        std::cout << "MtpSwitchFunc1Test (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpSwitchFunc1Test (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 3.2;
        rcut_smooth = 3.0;
    }

    void TearDown() override {
    }
};  // class : MtpSwitchFunc1Test


class MtpSwitchFunc2Test : public ::testing::Test
{
protected:
    double rcut;
    double rcut_smooth;
    double distance_ij;     // absolute value
    
    static void SetUpTestSuite() {
        std::cout << "MtpSwitchFunc2 (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpSwitchFunc2 (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 3.2;
        rcut_smooth = 3.0;
    }

    void TearDown() override {
    }
};  // class : SwitchFunc2Test


class ChebyshevPolyTest : public ::testing::Test {
protected:
    double rcut;
    double rcut_smooth;
    double distance_ij;
    static void SetUpTestSuite() {
        std::cout << "ChebyshevPolyTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "ChebyshevPolyTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 3.2;
        rcut_smooth = 3.0;
    }

    void TearDown() override {
    }
};  // class : ChebyshevPolyTest




TEST_F(MtpSwitchFunc1Test, copy_constructor) {
    distance_ij = 3.1;
    matersdk::mtp::MtpSwitchFunc1<double> msf1(rcut, rcut_smooth);
    matersdk::mtp::MtpSwitchFunc1<double> msf2(msf1);
}

TEST_F(MtpSwitchFunc1Test, assignment_operator) {
    distance_ij = 3.1;
    matersdk::mtp::MtpSwitchFunc1<double> msf1(rcut, rcut_smooth);
    matersdk::mtp::MtpSwitchFunc1<double> msf2(rcut-0.1, rcut_smooth);
    msf1 = msf2;
}

TEST_F(MtpSwitchFunc1Test, get_result_and_deriv2r) {
    distance_ij = 3.0;
    matersdk::mtp::MtpSwitchFunc1<double> msf1(rcut, rcut_smooth);
    ASSERT_DOUBLE_EQ(msf1.get_result(distance_ij), -1);
    printf("distance_ij = %5f, switch_func1_result = %6f, switch_func1_deriv2r = %6f\n",
        distance_ij,
        msf1.get_result(distance_ij),
        msf1.get_deriv2r());

    distance_ij = 3.1;
    matersdk::mtp::MtpSwitchFunc1<double> msf2(rcut, rcut_smooth);
    ASSERT_GE(msf2.get_result(distance_ij), -1);
    ASSERT_LE(msf2.get_result(distance_ij), 1);
    printf("distance_ij = %5f, switch_func1_result = %6f, switch_func1_deriv2r = %6f\n",
        distance_ij,
        msf2.get_result(distance_ij),
        msf2.get_deriv2r());

    distance_ij = 3.2;
    matersdk::mtp::MtpSwitchFunc1<double> msf3(rcut, rcut_smooth);
    ASSERT_DOUBLE_EQ(msf3.get_result(distance_ij), 1);
    printf("distance_ij = %5f, switch_func1_result = %6f, switch_func1_deriv2r = %6f\n",
        distance_ij,
        msf3.get_result(distance_ij),
        msf3.get_deriv2r());
}



TEST_F(MtpSwitchFunc2Test, copy_constructor) 
{
    distance_ij = 3.1;
    matersdk::mtp::MtpSwitchFunc2<double> msf1(rcut, rcut_smooth);
    matersdk::mtp::MtpSwitchFunc2<double> msf2(msf1);
}

TEST_F(MtpSwitchFunc2Test, assignment_operator)
{
    distance_ij = 3.1;
    matersdk::mtp::MtpSwitchFunc2<double> msf1(rcut, rcut_smooth);
    matersdk::mtp::MtpSwitchFunc2<double> msf2(rcut-0.1, rcut_smooth);
    msf1 = msf2;
}

TEST_F(MtpSwitchFunc2Test, get_result_and_deriv2r) {
    distance_ij = 3.0;
    matersdk::mtp::MtpSwitchFunc2<double> msf1(rcut, rcut_smooth);
    printf("distance_ij = %5f, switch_func2_result = %6f, switch_func1_deriv2r = %6f\n", 
        distance_ij, 
        msf1.get_result(distance_ij),
        msf1.get_deriv2r(distance_ij));
    ASSERT_DOUBLE_EQ(msf1.get_result(distance_ij), 1.0);

    distance_ij = 3.05;
    matersdk::mtp::MtpSwitchFunc2<double> msf2(rcut, rcut_smooth);
    printf("distance_ij = %5f, switch_func2_result = %6f, switch_func1_deriv2r = %6f\n", 
        distance_ij, 
        msf2.get_result(distance_ij),
        msf2.get_deriv2r(distance_ij));
    ASSERT_GE(msf2.get_result(distance_ij), 0.0);
    ASSERT_LE(msf2.get_result(distance_ij), 1.0);

    distance_ij = 3.2;
    matersdk::mtp::MtpSwitchFunc2<double> msf3(rcut, rcut_smooth);
    printf("distance_ij = %5f, switch_func2_result = %6f, switch_func1_deriv2r = %6f\n", 
        distance_ij, 
        msf3.get_result(distance_ij),
        msf3.get_deriv2r(distance_ij));
    ASSERT_DOUBLE_EQ(msf3.get_result(distance_ij), 0.0);
}


TEST_F(ChebyshevPolyTest, build) {
    distance_ij = 3.1;
    matersdk::mtp::ChebyshevPoly<double> chebyshev(8, rcut, rcut_smooth);
    chebyshev.build(distance_ij);
    for (int ii=0; ii<chebyshev.size(); ii++) {
        printf("%6f, ", chebyshev.get_result()[ii]);
    }
    printf("\n");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
