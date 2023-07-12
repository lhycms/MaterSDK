/**
 * @file CpuNeighborList.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2023-07-04
 * 
 * @copyright Copyright (c) 2023
 * 
 * @ref: 1. https://github.com/openmm/openmm/blob/master/platforms/cpu/include/CpuNeighborList.h
 *       2. https://github.com/openmm/openmm/blob/master/platforms/cpu/src/CpuNeighborList.cpp
 */
#ifndef CPU_NEIGHBOR_LIST_H
#define CPU_NEIGHBOR_LIST_H


#include <cmath>
#include <vector>

#include "../../../../core/include/vecx.h"
#include "../../../../core/include/hardware.h"
#include "../../../../core/include/AlignedArray.h"


namespace matersdk {

/**
 * @brief The Index of Voxel.
 * 
 */
class VoxelIndex {
public:
    VoxelIndex() : y(0), z(0)
    {};

    VoxelIndex(int y, int z): y(y), z(z)
    {};

    int y;
    int z;
};  // class VoxelIndex


class CpuNeighborList {
public:
    class Voxels {
    public:
        /**
         * @brief Construct a new Voxels object
         * 
         * @param blockSize 
         * @param vsy 
         * @param vsz 
         * @param miny 
         * @param maxy 
         * @param minz 
         * @param maxz 
         * @param boxVectors 
         * @param usePeriodic 
         */
        Voxels(int blockSize, 
            float vsy, float vsz, 
            float miny, float maxy, float minz, float maxz,
            const Vec3* boxVectors, bool usePeriodic);


        /**
         * @brief Get the Voxel Index object
         * 
         * @param location 
         * @return VoxelIndex 
         */
        VoxelIndex getVoxelIndex(const float* location) const;

    private:
        int blockSize;
        float voxelSizeY, voxelSizeZ;
        float miny, maxy, minz, maxz;
        int ny, nz;
        float periodicBoxSize[3], recipBoxSize[3];
        bool triclinic;
        float periodicBoxVectors[3][3];
        const bool usePeriodic;
        std::vector<std::vector<std::vector<std::pair<float, int>>>> bins;
    }; // class Voxels

    //class NeighborIterator;
    //CpuNeighborList(int blockSize);
    /**
     * @brief Compute the neighbor list based on the current positions of atoms
     * 
     * @param numAtoms          the number of atoms in the system
     * @param atomLocations     the positions of the atoms
     * 
     */


private:
    int blockSize;
    int num;
};



} /* namespace: matersdk */


#endif