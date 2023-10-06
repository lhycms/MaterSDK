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
    matersdk::mtp::Combinations combinations(mjus_njus_lst, false);
    combinations.show();
}


TEST_F(CombinationsTest, get_combinations) {
    matersdk::mtp::Combinations combinations(mjus_njus_lst, false);
    const std::vector<std::vector<std::pair<int, int>>> mjus_njus_lst = combinations.get_combinations();

    for (int ii=0; ii<mjus_njus_lst.size(); ii++) {
        printf("[mju, nju]: \t");
        int level = 0;
        for (int jj=0; jj<mjus_njus_lst[ii].size(); jj++) {
            printf(
                "[%4d, %4d],  ",
                this->mjus_njus_lst[ii][jj].first,
                this->mjus_njus_lst[ii][jj].second
            );
            level += matersdk::mtp::Combinations::get_level(this->mjus_njus_lst[ii][jj].first, this->mjus_njus_lst[ii][jj].second);
        }
        printf("level = %4d,\tno.%3d\n", level, ii);
    }
}


TEST_F(CombinationsTest, combinationsSortBasis_w1) {
    matersdk::mtp::Combinations combinations(mjus_njus_lst, false);
    matersdk::mtp::CombinationsSortBasis combination_sort_basis(combinations);

    printf("+++ %d\n", combination_sort_basis(15, 16));
}


TEST_F(CombinationsTest, combinationsSortBasis_w2) {
    matersdk::mtp::Combinations combinations(mjus_njus_lst, false);

    int* indices = (int*)malloc(sizeof(int) * combinations.get_num_combinations());
    for (int ii=0; ii<combinations.get_num_combinations(); ii++) {
        indices[ii] = ii;
    }
    std::sort(
            indices, 
            indices + combinations.get_num_combinations(), 
            matersdk::mtp::CombinationsSortBasis(combinations));
    for (int ii=0; ii<combinations.get_num_combinations(); ii++)
        printf("%4d, ", indices[ii]);
    printf("\n");
}


TEST_F(CombinationsTest, combinationsArrangement) {
    matersdk::mtp::Combinations combinations(mjus_njus_lst, false);

    // Step 1. get `int* new_indices`
    int* indices = (int*)malloc(sizeof(int) * combinations.get_num_combinations());
    for (int ii=0; ii<combinations.get_num_combinations(); ii++) {
        indices[ii] = ii;
    }
    std::sort(
            indices,
            indices + combinations.get_num_combinations(),
            matersdk::mtp::CombinationsSortBasis(combinations));
    
    // Step 2. Arrange `Combinations` according to `new_indices`
    matersdk::mtp::CombinationsArrangement combination_arrangement(combinations, indices);
    matersdk::mtp::Combinations sorted_combinations = combination_arrangement.arrange();

    // Step 3.
    sorted_combinations.show();
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}