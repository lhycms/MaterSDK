#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include "../include/mtpMAngularOp.h"


class MtpMAngularOpTest: public ::testing::Test {
protected:
    at::Tensor relative_coord_tensor;
    int nu_0;
    int nu_1;
    int nu_2;

    static void SetUpTestSuite() {
        std::cout << "MtpMAngularOpTest (TestSuite) is setting up...\n";        
    }

    static void TearDownTestSuite() {
        std::cout << "MtpMAngularOpTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        c10::TensorOptions options = c10::TensorOptions()
            .dtype(torch::kFloat32)
            .device(c10::kCPU);
        relative_coord_tensor = at::ones({3}, options);
        relative_coord_tensor.requires_grad_(true);
        nu_0 = 0;
        nu_1 = 1;
        nu_2 = 2;
    }

    void TearDown() override {
    }
};  // class : MtpMAngularOpTest


TEST_F(MtpMAngularOpTest, forward_and_backward) {
    at::Tensor mtp_angular_0 = matersdk::mtp::MtpMAngularOp(
        relative_coord_tensor,
        nu_0);
    at::Tensor mtp_angular_1 = matersdk::mtp::MtpMAngularOp(
        relative_coord_tensor,
        nu_1);
    at::Tensor mtp_angular_2 = matersdk::mtp::MtpMAngularOp(
        relative_coord_tensor,
        nu_2);
    assert( mtp_angular_0.dim() == 1 );
    assert( mtp_angular_1.dim() == 1 );
    assert( mtp_angular_2.dim() == 2 );
    //std::cout << mtp_angular_0.sizes() << std::endl;
    //std::cout << mtp_angular_1.sizes() << std::endl;
    //std::cout << mtp_angular_2.sizes() << std::endl;

    auto result = mtp_angular_2.sum();
    result.backward();
    assert (relative_coord_tensor.dim() == 1);
    //std::cout << relative_coord_tensor.grad() << std::endl;
}


TEST_F(MtpMAngularOpTest, speed) {
    int times = 5 * 1E5;
    auto time1 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int ii=0; ii<times; ii++) {
        at::Tensor mtp_angular = matersdk::mtp::MtpMAngularOp(
            relative_coord_tensor,
            2);
    }
    auto time2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1);
    std::cout << "Time costing in (ms) = " << duration.count() << std::endl;
}


TEST_F(MtpMAngularOpTest, deriv_accuracy) {
    at::Tensor mtp_angular_tensor = matersdk::mtp::MtpMAngularOp(
        relative_coord_tensor,
        2);
    auto result = mtp_angular_tensor.sum();
    result.backward();
    std::cout << "1. Partial derivative wrt. x calculated by Autograd = " << relative_coord_tensor.grad()[0].item<float>() << std::endl;

    float* relative_coord = relative_coord_tensor.data_ptr<float>();
    relative_coord[0] += 0.001;
    at::Tensor mtp_angular_tensor1 = matersdk::mtp::MtpMAngularOp(
        relative_coord_tensor,
        2);
    auto result1 = mtp_angular_tensor1.sum();
    std::cout << "2. Partial derivative wrt. x calculated by finite difference = " << (result1 - result).item<float>() / 0.001 << std::endl;
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}