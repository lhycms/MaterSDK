#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include "../include/mtpMAngularOp.h"


class MtpMAngularTest: public ::testing::Test {
protected:
    at::Tensor relative_coord_tensor;
    int nu_0;
    int nu_1;
    int nu_2;

    static void SetUpTestSuite() {
        std::cout << "MtpMAngularTest (TestSuite) is setting up...\n";        
    }

    static void TearDownTestSuite() {
        std::cout << "MtpMAngularTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        c10::TensorOptions options = c10::TensorOptions()
            .dtype(torch::kFloat32)
            .device(c10::kCPU);
        relative_coord_tensor = at::ones({3}, options);
        nu_0 = 0;
        nu_1 = 1;
        nu_2 = 2;
    }

    void TearDown() override {
    }
};  // class : MtpMAngularTest


TEST_F(MtpMAngularTest, all) {
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
}


TEST_F(MtpMAngularTest, speed) {
    int times = 1E5;
    auto time1 = std::chrono::high_resolution_clock::now();
    for (int ii=0; ii<times; ii++) {
        at::Tensor mtp_angular = matersdk::mtp::MtpMAngularOp(
            relative_coord_tensor,
            2);
    }
    auto time2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1);
    std::cout << "Time costing in (ms) = " << duration.count() << std::endl;
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}