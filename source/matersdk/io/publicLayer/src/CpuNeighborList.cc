#include "../include/CpuNeighborList.h"


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



class CpuNeighborList::Voxels {
public:
    /**
     * @brief 
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
            const Vec3* boxVectors, bool usePeriodic) :
        blockSize(blockSize), 
        voxelSizeY(vsy), voxelSizeZ(vsz),
        miny(miny), maxy(maxy), minz(minz), maxz(maxz),
        usePeriodic(usePeriodic) 
    {
        // Step 1. Initialize `this->periodicBoxVectors`
        for (int ii=0; ii<3; ii++) {
            for (int jj=0; jj<3; jj++) {
                float temp = (float)boxVectors[i][j];
                this->periodicBoxVectors[i][j] = temp;
            }
        }
        this->periodicBoxSize[0] = (float)boxVectors[0][0];
        this->periodicBoxSize[1] = (float)boxVectors[1][1];
        this->periodicBoxSize[2] = (float)boxVectors[2][2];

        // Step 2. Initialize `this->recipBoxVectors`
        this->recipBoxSize[0] = (float)(1.0 / boxVectors[0][0]);
        this->recipBoxSize[1] = (float)(1.0 / boxVectors[1][1]);
        this->recipBoxSize[2] = (float)(1.0 / boxVectors[2][2]);

        // Step 3. 判断Box是否是三斜的 (triclinic)
        triclinic = (
                boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
                boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
                boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0
            );

        // Step 4. 
        if (this->usePeriodic) {
            this->ny = (int) floorf(boxVectors[1][1] / this->voxelSizeY + 0.5f);
            this->nz = (int) floorf(boxVectors[2][2] / this->voxelSizeZ + 0.5f);
            this->voxelSizeY = boxVectors[1][1] / this->ny;
            this->voxelSizeY = boxVectors[2][2] / this->nz;
        } else {
            ny = NULL;
        }

        this->bins.resize(ny);
        for (int ii=0; ii<ny; ii++) {
            this->bins[ii].resize(nz);
        }

    } // class: Voxels


    /**
     * @brief Get the voxel index containing a particular location.
     * 
     * @param location 
     */
    VoxelIndex getVoxelIndex(const float* location) const {
        float yperiodic, zperiodic;
        if (!this->usePeriodic) {
            yperiodic = location[1] - miny;
            zperiodic = location[2] - minz;
        } else {
            
        }

    }



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
}; // class CpuNeighborList::Voxels



} // namespace: matersdk