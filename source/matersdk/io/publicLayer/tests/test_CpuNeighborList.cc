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
    matersdk::Vec3* boxVectors = new matersdk::Vec3[3];

    bool usePeriodic = true;


    static void SetUpTestSuite() {
        std::cout << "Start the CpuNeighborListTest (TestSuite)...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "End the CpuNeioghborListTest (TestSuite)...\n";
    }

    void SetUp() override {
        /**
         * @note You can just assign value in SetUp() function!
         * But not in protected region for pointer!
         * 
         */
        boxVectors[0][0] = 2.0;
        boxVectors[0][1] = 0.0;
        boxVectors[0][2] = 0.0;
        boxVectors[1][0] = 3.0;
        boxVectors[1][1] = 2.0;
        boxVectors[1][2] = 0.0;
        boxVectors[2][0] = 0.0;
        boxVectors[2][1] = 1.0;
        boxVectors[2][2] = 2.0;
        
    }

    void TearDown() override {
        delete [] boxVectors;
    }
};   // class CpuNeighborListTest



TEST_F(CpuNeighborListTest, getVoxelIndex) {
    matersdk::CpuNeighborList::Voxels voxels(
                                blockSize, 
                                vsy, vsz,
                                miny, maxy, minz, maxz,
                                boxVectors, usePeriodic);
    float location_1[3] = {1.0, 4.0, 12.0};
    matersdk::VoxelIndex index = voxels.getVoxelIndex(location_1);
    std::cout << index.y << ", " << index.z << std::endl;

    float location_2[3] = {1.0, 4.0, 11.0};
    matersdk::VoxelIndex index_2 = voxels.getVoxelIndex(location_2);
    std::cout << index_2.y << ", " << index_2.z << std::endl;
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}