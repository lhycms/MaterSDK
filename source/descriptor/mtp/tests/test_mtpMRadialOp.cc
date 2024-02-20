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
    at::Tensor rs_tensor;
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
        rs_tensor = at::zeros({3}, options);
        double* rs = rs_tensor.data_ptr<double>();
        rs[0] = rcut;
        rs[1] = rcut_smooth;
        rs[2] = distance_ij;
    }

    void TearDown() override {
    }
};  // class : MtpQOpTest


TEST_F(MtpQOpTest, apply) {
    at::Tensor result = matersdk::mtp::MtpQOp(size, rs_tensor)[0];
    std::cout << result << std::endl;
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
