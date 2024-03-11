#include <gtest/gtest.h>
#include <iostream>
#include <stdio.h>
#include "../include/mtpParam.h"


class MtpExceptionTest : public ::testing::Test
{
protected:
    static void SetUpTestSuite() {
        std::cout << "MtpException (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpException (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
    }

    void TearDown() override {
    }
};  // class : MtpExceptionTest

class MtpParamTest : public ::testing::Test
{
protected:
    std::string filename;

    static void SetUpTestSuite() {
        std::cout << "MtpParamTest (TestSuite) is setting up...\n";
    }

    static void TearDownTestSuite() {
        std::cout << "MtpParamTest (TestSuite) is tearing down...\n";
    }

    void SetUp() override {
        filename = "/data/home/liuhanyu/hyliu/code/matersdk/source/descriptor/mtpr/MTP_templates/06.almtp";
    }

    void TearDown() override {
    }
};  // class : MtpParam




TEST_F(MtpExceptionTest, throw_exception) {
    try {
        throw matersdk::mtpr::MtpException("gtest: throw mtp exception.");
    } catch (const matersdk::mtpr::MtpException& e) {
        std::cerr << e.what() << std::endl;
    }
}

TEST_F(MtpExceptionTest, MtpError) {
    using namespace matersdk::mtpr;

    try {
        MtpError("gtest: throw mtp exception.");
    } catch (const MtpException& e) {
        std::cerr << e.what() << std::endl;
    }
}

TEST_F(MtpParamTest, load) {
    matersdk::mtpr::MtpParam mtp_param;
    mtp_param._load(filename);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}