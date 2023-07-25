#ifndef CORE_VEC3_OPERATION_H
#define CORE_VEC3_OPERATION_H


#include <stdlib.h>
#include <cmath>

namespace matersdk {
namespace vec3Operation
{


template <typename CoordType>
CoordType* dot(CoordType* vec1, CoordType* vec2);

template <typename CoordType>
CoordType* cross(CoordType* vec1, CoordType* vec2);

template <typename CoordType>
CoordType* norm(CoordType* vec);




template <typename CoordType>
CoordType* cross(CoordType* vec1, CoordType* vec2) {
    CoordType* vertical_vec = (CoordType*)malloc(sizeof(CoordType) * 3);
    vertical_vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vertical_vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vertical_vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return vertical_vec;
}


template <typename CoordType>
CoordType* norm(CoordType* vec) {
    CoordType* unit_vec = (CoordType*)malloc(sizeof(CoordType) * 3);
    CoordType vec_length = std::sqrt( pow(vec[0], 2) + std::pow(vec[1], 2) + std::pow(vec[2], 2) );
    unit_vec[0] = vec[0] / vec_length;
    unit_vec[1] = vec[1] / vec_length;
    unit_vec[2] = vec[2] / vec_length;

    return unit_vec;
}


    
}   // namespace: vecOperation
}   // namespace: matersdk


#endif