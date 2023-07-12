#include <gtest/gtest.h>
#include <iostream>

// ./bin/matersdk/io/publicLayer/test_CpuNeighborList
#include "../include/CpuNeighborList.h"



class CpuNeighborListTest : public ::testing::Test {
protected:
    int blockSize = 10;
    float vsy = 1;
    float vsz = 1;
    float miny = 0;
    float maxy = 2;
    float minz = 0;
    float maxz = 2;
    const matersdk::Vec3* boxVectors = new matersdk::Vec3(1.0, 4.0, 12.0);
    bool usePeriodic = true;


    static void SetUpTestSuite() {
        std::cout << "Start the CpuNeighborListTest (TestSuite)...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "End the CpuNeioghborListTest (TestSuite)...\n";
    }

    void SetUp() override {

    }

    void TearDown() override {

    }
};   // class CpuNeighborListTest



TEST_F(CpuNeighborListTest, getVoxelIndex) {
    matersdk::CpuNeighborList::Voxels voxels(
                                blockSize, 
                                vsy, vsz,
                                miny, maxy, minz, maxz,
                                boxVectors, usePeriodic);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}