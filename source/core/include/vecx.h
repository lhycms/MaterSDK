// Copy from https://github.com/openmm/openmm/blob/5e9134005d3ca013572979289ce00ec32038e4f1/openmmapi/include/openmm/Vec3.h#L45
#ifndef CORE_VECX_H
#define CORE_VECX_H


#include <cassert>


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

    /**
     * @brief Create a Vec3 with specified x, y, and z components
     * 
     */
    Vec3(double x, double y, double z) {
        this->data[0] = x;
        this->data[1] = y;
        this->data[2] = z;
    }

    /**
     * @brief Take value using index
     * @param index int. 
     */
    double operator[](int index) const {
        assert(index>=0 && index<3);
        return this->data[index];
    }

    double& operator[](int index) {
        assert(index>=0 && index<3);
        return this->data[index];
    }

    /**
     * @brief Overload the `==` operator
     * 
     */
    bool operator==(const Vec3 &rhs) const {
        const Vec3 &lhs = *this;
        return (
                lhs[0] == rhs[0] && \
                lhs[1] == rhs[1] && \
                lhs[2] == rhs[2]
        );
    }

    /**
     * @brief Overload the `!=` operator 
     * 
     * @param rhs
     * @return 
     */
    bool operator!=(const Vec3 &rhs) const {
        const Vec3 &lhs = *this;
        return (
                lhs[0] != rhs[0] || \
                lhs[1] != rhs[1] || \
                lhs[2] != rhs[2]
        );
    }

    /**
     * @brief Overload the `+`(unary plus) operator
     * 
     */
    Vec3 operator+() const {
        const Vec3 &lhs = *this;
        return Vec3(lhs[0], lhs[1], lhs[2]);
    }


    /**
     * @brief Overload the `+`(binary plus) operator
     * 
     */
    Vec3 operator+(const Vec3 &rhs) const {
        const Vec3 &lhs = *this;
        return Vec3(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]);
    }

    /**
     * @brief Overload the `+=` operator
     * 
     */
    Vec3& operator+=(const Vec3 &rhs) {
        this->data[0] += rhs[0];
        this->data[1] += rhs[1];
        this->data[2] += rhs[2];
        return *this;
    }


private:
    double data[3];
};


}

#endif