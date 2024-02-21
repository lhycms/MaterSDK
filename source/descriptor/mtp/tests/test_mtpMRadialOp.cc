#include <gtest/gtest.h>
#include <stdio.h>
#include <iostream>
#include "../include/mtpMRadialOp.h"


class MtpQOpTest : public ::testing::Test
{
protected:
    int64_t size;
    double rcut;
    double rcut_smooth;
    double distance_ij;
    c10::TensorOptions options;
    at::Tensor rcuts_tensor;
    at::Tensor rcs_tensor;
    static void SetUpTestSuite() {
        std::cout << "MtpQOpTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpQOpTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        size = 8;
        rcut = 3.2;
        rcut_smooth = 3.0;
        distance_ij = 3.05;
        options = c10::TensorOptions()
            .dtype(torch::kFloat64)
            .device(c10::kCPU);
        rcuts_tensor = at::zeros({2}, options);
        double* rcuts = rcuts_tensor.data_ptr<double>();
        rcuts[0] = rcut;
        rcuts[1] = rcut_smooth;
        rcs_tensor = at::zeros({3, 3}, options);    // .shape = [nneighs, 3]
        double* rcs = rcs_tensor.data_ptr<double>();
        rcs[0*3 + 0] = 1.595158;    //3;
        rcs[0*3 + 1] = -0.920965;   //0;
        rcs[0*3 + 2] = -1.564884;   //0
        rcs[1*3 + 0] = 3.05;
        rcs[1*3 + 1] = 0;
        rcs[1*3 + 2] = 0;      
        rcs[2*3 + 0] = 0;           //3;
        rcs[2*3 + 1] = 0;
        rcs[2*3 + 2] = 0;
        rcs_tensor.requires_grad_(true);
    }

    void TearDown() override {
    }
};  // class : MtpQOpTest


TEST_F(MtpQOpTest, apply) {
    at::Tensor result = matersdk::mtp::MtpQOp(size, rcuts_tensor, rcs_tensor)[0];
    std::cout << "Result of Q(x) = \n" << result << std::endl;
    at::Tensor sum = result.sum();
    sum.backward();
    std::cout << "Partial derivative wrt. x_{ij} of second neigh =\n" << rcs_tensor.grad()[1][0] << std::endl;
    
    at::Tensor rcs_tensor_ = at::zeros({3, 3}, options);
    rcs_tensor_[0][0] = 3.0000;
    rcs_tensor_[0][1] = 0;
    rcs_tensor_[0][2] = 0;
    rcs_tensor_[1][0] = 3.0501;
    rcs_tensor_[1][1] = 0;
    rcs_tensor_[1][2] = 0;
    rcs_tensor_[2][0] = 3.2000;
    rcs_tensor_[2][1] = 0;
    rcs_tensor_[2][2] = 0;
    rcs_tensor_.requires_grad_(true);
    at::Tensor result_ = matersdk::mtp::MtpQOp(size, rcuts_tensor, rcs_tensor_)[0];
    at::Tensor sum_ = result_.sum();
    std::cout << "Result of Q(x+/delta{x}}) = \n" << result_ << std::endl;
    std::cout << "Partial derivative wrt. x_{ij} of second neigh =\n" << ((sum_ - sum) / 0.0001) << std::endl;
    
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
