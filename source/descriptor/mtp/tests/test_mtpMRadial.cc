#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include "../include/mtpMRadial.h"


class SwitchFunc1Test : public ::testing::Test
{
protected:
    double rcut;
    double rcut_smooth;
    double distance_ij;    // absolute value
    static void SetUpTestSuite() {
        std::cout << "SwitchFunc1Test (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "SwitchFunc1Test (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rcut = 3.2;
        rcut_smooth = 3.0;
    }

    void TearDown() override {
    }
};  // class : SwitchFunc1Test


TEST_F(SwitchFunc1Test, get_result) {
    distance_ij = 3.0;
    matersdk::mtp::MtpSwitchFunc1<double> msf1(rcut, rcut_smooth);
    ASSERT_DOUBLE_EQ(msf1.get_result(distance_ij), -1);

    distance_ij = 3.2;
    matersdk::mtp::MtpSwitchFunc1<double> msf3(rcut, rcut_smooth);
    ASSERT_DOUBLE_EQ(msf3.get_result(distance_ij), 1);
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
