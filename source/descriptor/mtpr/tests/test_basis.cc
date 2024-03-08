#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include "../include/basis.h"


class SwitchFunctionTest : public ::testing::Test
{
protected:
    double rmax;
    double rmin;
    double distance_ij;

    static void SetUpTestSuite() {
        std::cout << "SwitchFunctionTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSutie() {
        std::cout << "SwitchFunctionTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        rmax = 5.0;
        rmin = 2.0;
        distance_ij = 3.14;
    }

    void TearDown() override {
    }
};  // class : SwitchFunctionTest


class RB_ChebyshevTest : public ::testing::Test
{
protected:
    int size;
    double rmax;
    double rmin;
    double distance_ij;

    static void SetUpTestSuite() {
        std::cout << "RB_ChebyshevTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "RB_ChebyshevTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        size = 8;
        rmax = 5.0;
        rmin = 2.0;
    }

    void TearDown() override {
    }
};  // class : RB_ChebyshevTest




TEST_F(SwitchFunctionTest, init)
{
    matersdk::mtpr::SwitchFunction<double> swf(rmax, rmin);
    distance_ij = rmin;
    ASSERT_DOUBLE_EQ(swf.value(distance_ij), 1);
    ASSERT_DOUBLE_EQ(swf.der2r(distance_ij), 0);

    distance_ij = rmax;
    ASSERT_DOUBLE_EQ(swf.value(distance_ij), 0);
    ASSERT_DOUBLE_EQ(swf.der2r(distance_ij), 0);
}

TEST_F(SwitchFunctionTest, der_accuracy)
{
    distance_ij = 3.14;
    matersdk::mtpr::SwitchFunction<double> swf(rmax, rmin);

    double der2r = swf.der2r(distance_ij);
    double value1 = swf.value(distance_ij);
    double value2 = swf.value(distance_ij + 0.0001);
    double der2r_ = (value2 - value1) / 0.0001;

std::cout << "Custom code method: Deriv wrt. r = " << der2r << std::endl;
std::cout << "Finite difference method: Deriv wrt. r = " << der2r_ << std::endl;
}


TEST_F(RB_ChebyshevTest, build) 
{
    matersdk::mtpr::RB_Chebyshev<double>* rb_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    ASSERT_EQ(rb_ptr->size(), size);
    ASSERT_EQ(rb_ptr->rmax(), rmax);
    ASSERT_EQ(rb_ptr->rmin(), rmin);

    distance_ij = 3.14;
    rb_ptr->build(distance_ij);
//rb_ptr->show();
    delete rb_ptr;
}

TEST_F(RB_ChebyshevTest, der_accuracy)
{
    distance_ij = 3.14;
    matersdk::mtpr::RB_Chebyshev<double>* rb1_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    rb1_ptr->build(distance_ij);

    double distance_ij_ = 3.14 + 0.0001;
    matersdk::mtpr::RB_Chebyshev<double>* rb2_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    rb2_ptr->build(distance_ij_);

printf("Custom code method: deriv of RB_Chebyshev wrt. r:\n\t");
for (int ii=0; ii<rb1_ptr->size(); ii++)
    printf("%10lf, ", rb1_ptr->ders2r()[ii]);
printf("\n");
printf("Finite difference method: deriv of RB_Chebyshev wrt. r:\n\t");
for (int ii=0; ii<rb2_ptr->size(); ii++) {
    double tmp_der2r = (rb2_ptr->vals()[ii] - rb1_ptr->vals()[ii]) / 0.0001;
    printf("%10lf, ", tmp_der2r);
}
printf("\n");

    delete rb1_ptr;
    delete rb2_ptr;
}

TEST_F(RB_ChebyshevTest, copy_constructor)
{
    matersdk::mtpr::RB_Chebyshev<double>* rb1_ptr  = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    matersdk::mtpr::RB_Chebyshev<double> rb2(*rb1_ptr);

    ASSERT_EQ(rb1_ptr->size(), rb2.size());
    ASSERT_DOUBLE_EQ(rb1_ptr->rmin(), rb2.rmin());
    ASSERT_DOUBLE_EQ(rb1_ptr->rmax(), rb2.rmax());
    for (int ii=0; ii<rb1_ptr->size(); ii++) {
        ASSERT_DOUBLE_EQ(rb1_ptr->vals()[ii], rb2.vals()[ii]);
        ASSERT_DOUBLE_EQ(rb1_ptr->ders2uu()[ii], rb2.ders2uu()[ii]);
        ASSERT_DOUBLE_EQ(rb1_ptr->ders2r()[ii], rb2.ders2r()[ii]);
    }

    delete rb1_ptr;
}

TEST_F(RB_ChebyshevTest, copy_constructor_move)
{
    matersdk::mtpr::RB_Chebyshev<double>* rb1_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    matersdk::mtpr::RB_Chebyshev<double> rb2(std::move(*rb1_ptr));
    
    ASSERT_EQ(rb1_ptr->vals(), nullptr);
    ASSERT_EQ(rb1_ptr->ders2uu(), nullptr);
    ASSERT_EQ(rb1_ptr->ders2r(), nullptr);

    delete rb1_ptr;
}

TEST_F(RB_ChebyshevTest, assignment_operator)
{
    matersdk::mtpr::RB_Chebyshev<double>* rb1_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    matersdk::mtpr::RB_Chebyshev<double>* rb2_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin+1);

    *rb2_ptr = *rb1_ptr;

    ASSERT_EQ(rb1_ptr->size(), rb2_ptr->size());
    ASSERT_DOUBLE_EQ(rb1_ptr->rmax(), rb2_ptr->rmax());
    ASSERT_DOUBLE_EQ(rb1_ptr->rmin(), rb2_ptr->rmin());
    for (int ii=0; ii<rb1_ptr->size(); ii++) {
        ASSERT_DOUBLE_EQ(rb1_ptr->vals()[ii], rb2_ptr->vals()[ii]);
        ASSERT_DOUBLE_EQ(rb1_ptr->ders2uu()[ii], rb2_ptr->ders2uu()[ii]);
        ASSERT_DOUBLE_EQ(rb1_ptr->ders2r()[ii], rb2_ptr->ders2r()[ii]);
    }

    delete rb1_ptr;
    delete rb2_ptr;
}

TEST_F(RB_ChebyshevTest, assignment_operator_move)
{
    matersdk::mtpr::RB_Chebyshev<double>* rb1_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin);
    matersdk::mtpr::RB_Chebyshev<double>* rb2_ptr = new matersdk::mtpr::RB_Chebyshev<double>(size, rmax, rmin+1);

    *rb2_ptr = std::move(*rb1_ptr);

    ASSERT_EQ(rb1_ptr->vals(), nullptr);
    ASSERT_EQ(rb1_ptr->ders2uu(), nullptr);
    ASSERT_EQ(rb1_ptr->ders2r(), nullptr);
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}