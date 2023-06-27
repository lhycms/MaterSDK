#ifndef NEIGH_H
#define NEIGH_H


#include <pair>
#include <vector>


typedef std::pair<const Vec3*, AtomIndex> VoxelItem;
typedef std::vector<VoxelItem> Voxel;


// Part 1. NeighborListV1 declaration
class NeighborListV1 {
public:
    // Step 1.1. Constructor
    NeighborList();
}

#endif