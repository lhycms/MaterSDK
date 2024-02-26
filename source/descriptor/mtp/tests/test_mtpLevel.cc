#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include <utility>
#include "../include/mtpLevel.h"


class MtpMCoeffPairTest : public ::testing::Test 
{
protected:
    int mu_0, nu_0;
    std::pair<int, int> coeff_pair_0;

    static void SetUpTestSuite() {
        std::cout << "MtpMCoeffPairTest is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpMCoeffPairTest is tearing down...\n";
    }

    void SetUp() override {
        mu_0 = 1;
        nu_0 = 2;
        coeff_pair_0 = std::pair<int, int>(mu_0, nu_0);
    }

    void TearDown() override {
    }
};  // class : MtpMCoeffPairTest



TEST_F(MtpMCoeffPairTest, constructor_default)
{
    matersdk::mtp::MtpMCoeffPair mcp;
    ASSERT_EQ(mcp.level(), 0);
}

TEST_F(MtpMCoeffPairTest, constructor_1)
{
    matersdk::mtp::MtpMCoeffPair mcp(mu_0, nu_0);
    //std::cout << mcp.level() << std::endl;
    ASSERT_EQ(mcp.level(), 8);
}

TEST_F(MtpMCoeffPairTest, constructor_2)
{
    matersdk::mtp::MtpMCoeffPair mcp(coeff_pair_0);
    //std::cout << mcp.level() << std::endl;
    ASSERT_EQ(mcp.level(), 8);
}

TEST_F(MtpMCoeffPairTest, copy_constructor)
{
    matersdk::mtp::MtpMCoeffPair mcp_1(coeff_pair_0);
    matersdk::mtp::MtpMCoeffPair mcp_2(mcp_1);
    ASSERT_EQ(mcp_1.level(), mcp_2.level());
}

TEST_F(MtpMCoeffPairTest, assignment_operator)
{
    matersdk::mtp::MtpMCoeffPair mcp_1(coeff_pair_0);
    matersdk::mtp::MtpMCoeffPair mcp_2;
    ASSERT_EQ(mcp_1.level(), 8);
    ASSERT_EQ(mcp_2.level(), 0);
    mcp_2 = mcp_1;
    ASSERT_EQ(mcp_1.level(), 8);
    ASSERT_EQ(mcp_2.level(), 8);
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}