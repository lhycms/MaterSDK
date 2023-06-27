// Copy from https://github.com/openmm/openmm/blob/5e9134005d3ca013572979289ce00ec32038e4f1/openmmapi/include/openmm/Vec3.h#L45
#ifndef VECX_H
#define VECX_H

#include <cassert>


class Vec3 {
public:
    /*
        Description
        -----------
            1. Create a Vec3 whose elements are all 0.
    */
    Vec3() {
        data[0] = data[1] = data[2] = 0.0;
    }


    /*
        Description
        -----------
            1.Create a Vec3 with specified x, y and z components.
    */
    Vec3(double x, double y, double z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }


    /*
        Description
        -----------
            1. 
    */
    double operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    /*
        Description
        -----------
            1. 
    */
    double operator[](int index) {
        assert (index >= 0 && index < 3);
        return data[index];
    }

    /*
        Description
        -----------
            1. 
    */
    bool operator==(const Vec3& rhs) const {
        return (
            (this->data[0] == rhs[0]) && 
            (this->data[1] == rhs[1]) &&
            (this->data[2] == rhs[2])
        );
    }

    /*
        Description
        -----------
            1. 
    */
    bool operator!=(const Vec3& rhs) const {
        return (
            (this->data[0] != rhs[0]) ||
            (this->data[1] != rhs[1]) ||
            (this->data[2] != rhs[2])
        );
    }

    /*
        Description
        -----------
            1. Unary plus
    */
    Vec3 operator+() const {
        return Vec3(*this);
    }

    /*
        Description
        -----------
            1. plus
    */
    Vec3 operator+(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }

    /*
        Description
        -----------
            1. 
    */
    Vec3& operator+=(const Vec3& rhs) {
        this->data[0] += rhs[0];
        this->data[1] += rhs[1];
        this->data[2] += rhs[2];
        return *this;
    }

    /*
        Description
        -----------
            1. Unary minus
    */
    Vec3 operator-() const {
        const Vec3& lhs = *this;
        return Vec3(-lhs[0], -lhs[1], -lhs[2]);
    }

    /*
        Description
        -----------
            1. 
    */
    



private:
    double data[3];
}


#endif