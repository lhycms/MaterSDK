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
#include <pair>

#include "../../../../core/include/vecx.h"
#include "../../../../core/include/hardware.h"
#include "../../../../core/include/AlignedArray.h"


namespace matersdk {

class CpuNeighborList {
public:
    class Voxels;
    class NeighborIterator;
    CpuNeighborList(int blockSize);
    /**
     * @brief Compute the neighbor list based on the current positions of atoms
     * 
     * @param numAtoms          the number of atoms in the system
     * @param atomLocations     the positions of the atoms
     * 
     */


private:
    int blockSize;
};



} /* namespace: matersdk */


#endif