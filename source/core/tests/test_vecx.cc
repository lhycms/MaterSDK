#include <gtest/gtest.h>
#include <iostream>

// cmake -DBUILD_TEST=1 ..; make -j 8; ./bin/core/test_vecx 
#include "../include/vecx.h"


class Vec3Test : public ::testing::Test {
protected:
    matersdk::Vec3 *ptr_vec3_1;
    matersdk::Vec3 *ptr_vec3_2;

    static void SetUpTestSuite() {
        std::cout << "Set up Vec3Test (TestSuite)...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "Tear down Vec3Test (TestSuite)...\n";
    }

    void SetUp() override {
        ptr_vec3_1 = new matersdk::Vec3(1.0, 2.0, 3.0);
        ptr_vec3_2 = new matersdk::Vec3(1.0, 2.0, 3.0);
    }

    void TearDown() override {
        delete ptr_vec3_1;
        delete ptr_vec3_2;
    }
};


/**
 * @brief Construct a new test f object 
 * for `matersdk::Vec3::operator[]`
 * 
 */
TEST_F(Vec3Test, Index) {
    EXPECT_EQ((*ptr_vec3_1)[0], 1.0);
    EXPECT_EQ((*ptr_vec3_1)[1], 2.0);
    EXPECT_EQ((*ptr_vec3_1)[2], 3.0);


    //const matersdk::Vec3 vec3_const(1.0, 2.0, 3.0);
    //double &result = vec3_const[0];
    //result += 1;
    //EXPECT_EQ(vec3_const[0], 1.0);

    matersdk::Vec3 vec3_nonconst(1.0, 2.0, 3.0);
    double &result = vec3_nonconst[0];
    result += 1;
    EXPECT_EQ(vec3_nonconst[0], 2.0);
}

/**
 * @brief Construct a new test f object
 * for `matersdk::Vec3::operator==` and `matersdk::Vec3::operator!=`
 * 
 */
TEST_F(Vec3Test, EqNe) {
    EXPECT_EQ(*(ptr_vec3_1), *(ptr_vec3_2));

    matersdk::Vec3 vec3_3(2.0, 3.0, 4.0);
    EXPECT_NE(*(ptr_vec3_1), vec3_3);
}


/**
 * @brief Construct a new test f object
 * for `matersdk::Vec3::operator+` (Unary plus)
 * 
 */
TEST_F(Vec3Test, UnaryPlus) {
    EXPECT_EQ(+(*ptr_vec3_1), (*ptr_vec3_1));
}


/**
 * @brief Construct a new test f object
 * for `matersdk::Vec3::operator+` (Binary plus)
 * 
 */
TEST_F(Vec3Test, BinaryPlus) {
    matersdk::Vec3 result = (*ptr_vec3_1) + (*ptr_vec3_2);
    EXPECT_EQ(result, matersdk::Vec3(2.0, 4.0, 6.0));
}

/**
 * @brief Construct a new test f object
 * for `matersdk::Vec3::operator+=`
 * 
 */
TEST_F(Vec3Test, SelfPlus) {
    (*ptr_vec3_1) += (*ptr_vec3_1);
    EXPECT_EQ((*ptr_vec3_1), matersdk::Vec3(2.0, 4.0, 6.0));
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}