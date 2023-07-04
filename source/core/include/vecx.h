// Copy from https://github.com/openmm/openmm/blob/5e9134005d3ca013572979289ce00ec32038e4f1/openmmapi/include/openmm/Vec3.h#L45
#ifndef CORE_VECX_H
#define CORE_VECX_H

namespace matersdk {

/**
 * @brief This Class represents a three component vector. It is used for storing positions, 
 * velocities, and forces 
 */
class Vec3 {
public:
    /**
     * @brief Create a Vec3 whose elements are all 0.
     * 
     */
    Vec3() {
        this->data[0] = this->data[1] = this->data[2] = 0.0;
    }


private:
    double data[3];
}


}

#endif