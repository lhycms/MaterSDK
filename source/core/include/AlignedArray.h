/**
 * @file AlignedArray.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2023-07-04
 * 
 * @copyright Copyright (c) 2023
 * 
 * @ref 1. https://github.com/openmm/openmm/blob/644dc1ecc9e95b3c8e831803bb3e2ac925999f74/platforms/cpu/include/AlignedArray.h#L4
 */
#ifndef CORE_ALIGNED_ARRAY_H
#define CORE_ALIGNED_ARRAY_H

namespace matersdk {

/**
 * @brief 
 * 1. This class represents an array in memory whose `starting point is guaranted to
 * be aligned with a 16 bytes boundary` -- `it means that the memory address where the array begins 
 * is divisible by 16`. This can improve the performance of vectorized code, since loads and 
 * stores are more efficient.
 * 2. 通过将数组与 16 字节边界对齐，可以确保在执行向量化操作时，数据可以有效地加载到 SIMD 寄存器中。
 * 3. SIMD 指令设计用于同时操作多个数据元素，并且它们通常要求数据与特定边界对齐以获得最佳性能。
 * 4. In th case of `16-byte boundary`, it means that memory address divisible by 16 are 
 * considered aligned address.
 */
template <typename T>
class AlignedArray {
public:
    /**
     * @brief Default constructor, to allow AlignedArrays to be used inside collections
     * 
     */
    AlignedArray() : dataSize(0), baseData(0), data(0)
    {};

    /**
     * @brief Create an Aligned array that contains a specified number of elements.
     * 
     */
    AlignedArray(int size) {
        this->allocate(size);
    }

    /**
     * @brief Destructor
     * 
     */
    ~AlignedArray() {
        if (baseData != 0) {
            this->dataSize = 0;
            delete[] this->baseData;
        }
    }

    /**
     * @brief Get the number of elements in the array.
     */
    int size() const {
        return this->dataSize;
    }

private:
    /**
     * @brief This function is responsible for allocating memory for the aligned array.
     * 
     * @param size representing the number of elements to allocate
     */
    void allocate(int size) {
        this->dataSize = size;  // `this->dataSize` member variable stores the provided size.
        /* 
            * The size of allocated memory is calculated as `size*sizeof(T)+16`.
            * The additional 16 bytes ensure that starting address of the array will be aligned to 
            * a 16-byte boundary.
        */
        this->baseData = new char[size*sizeof(T)+16];
        /*
            * `offsetData` is set to point to the allocated memory plus 15 bytes.
            * This is done to ensure that the offset is `a multiple of 16`.
        */
        char *offsetData = this->baseData + 15;
        /*
            * The expression `(long long)offsetData & 0xF` calculates the lower 4 bits
        */
        offsetData -= (long long)offsetData & 0xF;
        this->data = (T*) offsetData;
    }
    int dataSize;
    char* baseData;
    T* data;

}; /* class: AlignedArray */


} /* namespace: matersdk */


#endif