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
    at::Tensor distances_tensor;
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
        distances_tensor = at::zeros({3}, options);
        double* distances = distances_tensor.data_ptr<double>();
        distances[0] = 3;
        distances[1] = 3.05;
        distances[2] = 3.2;
        distances_tensor.requires_grad_(true);
    }

    void TearDown() override {
    }
};  // class : MtpQOpTest


TEST_F(MtpQOpTest, apply) {
    at::Tensor result = matersdk::mtp::MtpQOp(size, rcuts_tensor, distances_tensor)[0];
    std::cout << result << std::endl;
    at::Tensor sum = result.sum();
    sum.backward();
    std::cout << distances_tensor.grad() << std::endl;
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
