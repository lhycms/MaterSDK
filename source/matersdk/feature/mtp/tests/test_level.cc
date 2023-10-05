#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include "../include/level.h"


class CombinationsTest: public ::testing::Test {
protected:
    std::vector<std::vector<std::pair<int, int>>> mjus_njus_lst;

    static void SetUpTestSuite() {
        std::cout << "CombinationsTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "CombinationsTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        mjus_njus_lst.resize(17);
        for (int ii=0; ii<mjus_njus_lst.size(); ii++) {
            mjus_njus_lst[ii].resize(2);
        }

        mjus_njus_lst[0][0].first = 0;
        mjus_njus_lst[0][0].second = 0;
        mjus_njus_lst[0][1].first = 0;
        mjus_njus_lst[0][1].second = 0;

        mjus_njus_lst[1][0].first = 0;
        mjus_njus_lst[1][0].second = 0;
        mjus_njus_lst[1][1].first = 0;
        mjus_njus_lst[1][1].second = 1;

        mjus_njus_lst[2][0].first = 0;
        mjus_njus_lst[2][0].second = 0;
        mjus_njus_lst[2][1].first = 0;
        mjus_njus_lst[2][1].second = 2;

        mjus_njus_lst[3][0].first = 0;
        mjus_njus_lst[3][0].second = 0;
        mjus_njus_lst[3][1].first = 0;
        mjus_njus_lst[3][1].second = 3;

        mjus_njus_lst[4][0].first = 0;
        mjus_njus_lst[4][0].second = 0;
        mjus_njus_lst[4][1].first = 0;
        mjus_njus_lst[4][1].second = 4;

        mjus_njus_lst[5][0].first = 0;
        mjus_njus_lst[5][0].second = 0;
        mjus_njus_lst[5][1].first = 1;
        mjus_njus_lst[5][1].second = 0;

        mjus_njus_lst[6][0].first = 0;
        mjus_njus_lst[6][0].second = 1;
        mjus_njus_lst[6][1].first = 0;
        mjus_njus_lst[6][1].second = 0;

        mjus_njus_lst[7][0].first = 0;
        mjus_njus_lst[7][0].second = 1;
        mjus_njus_lst[7][1].first = 0;
        mjus_njus_lst[7][1].second = 1;

        mjus_njus_lst[8][0].first = 0;
        mjus_njus_lst[8][0].second = 1;
        mjus_njus_lst[8][1].first = 0;
        mjus_njus_lst[8][1].second = 2;

        mjus_njus_lst[9][0].first = 0;
        mjus_njus_lst[9][0].second = 1;
        mjus_njus_lst[9][1].first = 0;
        mjus_njus_lst[9][1].second = 3;

        mjus_njus_lst[10][0].first = 0;
        mjus_njus_lst[10][0].second = 2;
        mjus_njus_lst[10][1].first = 0;
        mjus_njus_lst[10][1].second = 0;

        mjus_njus_lst[11][0].first = 0;
        mjus_njus_lst[11][0].second = 2;
        mjus_njus_lst[11][1].first = 0;
        mjus_njus_lst[11][1].second = 1;

        mjus_njus_lst[12][0].first = 0;
        mjus_njus_lst[12][0].second = 2;
        mjus_njus_lst[12][1].first = 0;
        mjus_njus_lst[12][1].second = 2;

        mjus_njus_lst[13][0].first = 0;
        mjus_njus_lst[13][0].second = 3;
        mjus_njus_lst[13][1].first = 0;
        mjus_njus_lst[13][1].second = 0;

        mjus_njus_lst[14][0].first = 0;
        mjus_njus_lst[14][0].second = 3;
        mjus_njus_lst[14][1].first = 0;
        mjus_njus_lst[14][1].second = 1;

        mjus_njus_lst[15][0].first = 0;
        mjus_njus_lst[15][0].second = 4;
        mjus_njus_lst[15][1].first = 0;
        mjus_njus_lst[15][1].second = 0;

        mjus_njus_lst[16][0].first = 1;
        mjus_njus_lst[16][0].second = 0;
        mjus_njus_lst[16][1].first = 0;
        mjus_njus_lst[16][1].second = 0;
    }

    void TearDown() override {

    }

};  // class : Combinations



TEST_F(CombinationsTest, constructor_1) {
    Combinations combinations(mjus_njus_lst, true);
    combinations.show();
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}