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
        std::cout << "MtpMCoeffPairTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpMCoeffPairTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        mu_0 = 1;
        nu_0 = 2;
        coeff_pair_0 = std::pair<int, int>(mu_0, nu_0);
    }

    void TearDown() override {
    }
};  // class : MtpMCoeffPairTest


class MtpMCoeffPairCombsTest : public ::testing::Test
{
protected:
    int max_level_1;
    int max_level_2;
    int aim_level_1;
    int aim_level_2;

    static void SetUpTestSuite() {
        std::cout << "MtpMCoeffPairCombs (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpMCoeffPairCombs (TestSuite) is tearding down...\n";
    }

    void SetUp() override {
        max_level_1 = 6;
        max_level_2 = 8;
        aim_level_1 = 6;
        aim_level_2 = 8;
    }

    void TearDown() override {
    }
};  // class : MtpMCoeffPairCombsTest


std::ostream& operator<<(
    std::ostream& COUT, 
    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> coeff_pair_combs)
{
    int count = 0;
    for (auto& coeff_pair_comb : coeff_pair_combs) {
        printf("Comb#%5d:\n\t", count);
        for (auto& coeff_pair : coeff_pair_comb) {
            printf("[%3d, %3d], ", coeff_pair.coeff_pair().first, coeff_pair.coeff_pair().second);
        }
        printf("\n");
        count++;
    }
    printf("MtpM Level = %3d, %5d combinations in total.", coeff_pair_combs.size(), coeff_pair_combs.size());
    return COUT;
}



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


TEST_F(MtpMCoeffPairCombsTest, constructor_default)
{
    matersdk::mtp::MtpMCoeffPairCombs mcpcs;
    //std::cout << mcpcs.coeff_pair_combs().size() << std::endl;
    ASSERT_EQ(mcpcs.coeff_pair_combs().size(), 0);
}

TEST_F(MtpMCoeffPairCombsTest, constructor_1)
{
    matersdk::mtp::MtpMCoeffPairCombs mcpcs_1(aim_level_1);
    //mcpcs_1.show();
    ASSERT_EQ(mcpcs_1.coeff_pair_combs().size(), 5);

    matersdk::mtp::MtpMCoeffPairCombs mcpcs_2(aim_level_2);
    mcpcs_2.show();
    ASSERT_EQ(mcpcs_2.coeff_pair_combs().size(), 9);
}

TEST_F(MtpMCoeffPairCombsTest, get_all_schemes4lev)
{
    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> combs_0 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_all_schemes4lev(0, 0, 0);
    ASSERT_EQ(combs_0.size(), 1);   // Contains empty scheme.

    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> combs_1 =
        matersdk::mtp::MtpMCoeffPairCombs::get_all_schemes4lev(aim_level_1, 0, 0);
    //std::cout << combs_1 << std::endl;
    ASSERT_EQ(combs_1.size(), 5);

    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> combs_2 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_all_schemes4lev(aim_level_2, 0, 0);
//std::cout << combs_2 << std::endl;
    ASSERT_EQ(combs_2.size(), 9);
}

TEST_F(MtpMCoeffPairCombsTest, get_contracted_combs)
{
    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> coeff_pair_combs_6 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_all_schemes4lev(aim_level_1, 0, 0);
    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> contracted_coeff_pair_combs_6 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_contracted_combs(coeff_pair_combs_6);
    //std::cout << contracted_coeff_pair_combs_6 << std::endl;
    ASSERT_EQ(contracted_coeff_pair_combs_6.size(), 3);

    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> coeff_pair_combs_8 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_all_schemes4lev(aim_level_2, 0, 0);
    std::vector<std::vector<matersdk::mtp::MtpMCoeffPair>> contracted_coeff_pair_combs_8 = 
        matersdk::mtp::MtpMCoeffPairCombs::get_contracted_combs(coeff_pair_combs_8);
    //std::cout << contracted_coeff_pair_combs_8 << std::endl;
    ASSERT_EQ(contracted_coeff_pair_combs_8.size(), 4);
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}