// Copy from https://github.com/openmm/openmm/blob/5e9134005d3ca013572979289ce00ec32038e4f1/openmmapi/include/openmm/Vec3.h#L45
#ifndef VECX_H
#define VECX_H

#include <iostream>
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
            1. Create a Vec3 with specified x, y and z components.
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
            1. Minus
    */
    Vec3 operator-(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
    }

    /*
        Description
        -----------
            1. 
    */
    Vec3& operator-=(const Vec3& rhs) {
        this->data[0] -= rhs[0];
        this->data[1] -= rhs[1];
        this->data[2] -= rhs[2];
    }

    /*
        Description
        -----------
            1. Scalar product
    */
    Vec3 operator*(double rhs) const {
        const Vec3& lhs = *this;
        return Vec3(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
    }

    /*
        Description
        -----------
            1. 
    */
    Vec3& operator*=(double rhs) {
        this->data[0] *= rhs;
        this->data[1] *= rhs;
        this->data[2] *= rhs;
        return *this;
    }
    
    /*
        Description
        -----------
            1. Scalar division
    */
    Vec3 operator/(double rhs) const {
        const Vec3& lhs = *this;
        double scale_factor = 1.0 / rhs;
        return Vec3(lhs[0]*scale_factor, lhs[1]*scale_factor, lhs[2]*scale_factor);
    }

    /*
        Description
        -----------
            1.
    */
    Vec3& operator/(double rhs) {
        double scale_factor = 1.0 / rhs;
        this->data[0] *= scale_factor;
        this->data[1] *= scale_factor;
        this->data[2] *= scale_factor;
        return *this;
    }

    /*
        Description
        -----------
            1. Dot production
    */
    double dot(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
    }

    /*
        Description
        -----------
            1. Cross production
    */
    Vec3 cross(const Vec3& rhs) const {
        const Vec3& lhs = *this;
        return Vec3(
                lhs[1]*rhs[2] - lhs[2]*rhs[1],
                lhs[2]*rhs[0] - lhs[0]*rhs[2],
                lhs[0]*rhs[1] - lhs[1]*rhs[0]
        );
    }


private:
    double data[3];
}


Vec3 operator*(double lhs, Vec3 rhs) {
    return Vec3(rhs[0]*lhs, rhs[1]*lhs, rhs[2]*lhs);
}


std::ostream& operator<<(std::ostream& COUT, const Vec3& vec3) {
    COUT << "[" << vec3[0] << ", " << vec3[1] << ", " << vec3[2] << "]";
    return COUT;
}


#endif