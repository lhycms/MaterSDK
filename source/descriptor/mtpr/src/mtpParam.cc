#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cstdlib>
#include "../include/mtpParam.h"


namespace matersdk {
namespace mtpr {

MtpParam::MtpParam()
{
    this->_alpha_moments_count = 0;
    this->_alpha_index_basic_count = 0;
    this->_alpha_index_basic = nullptr;
    this->_alpha_index_times_count = 0;
    this->_alpha_index_times = nullptr;
    this->_alpha_scalar_moments = 0;
    this->_alpha_moment_mapping = nullptr;
}

MtpParam::MtpParam(const std::string& filename)
{
    this->_load(filename);
}

void MtpParam::_load(const std::string& filename)
{
    std::ifstream ifs;
    ifs.open(filename);
    if ( !ifs.is_open() ) {
        MtpError((std::string)"Cannot open file " + filename);
    }
    char tmpline[1000];
    std::string tmpstr;

    for (int ii=0; ii<10; ii++) 
        ifs.ignore(5000, '\n');
    
    // alpha_moments_count
    ifs >> tmpstr;
    if (tmpstr != "alpha_moments_count")
        MtpError("Cannot read alpha_moment_counts.");
    ifs.ignore(2);
    ifs >> this->_alpha_moments_count;
    if (ifs.fail()) 
        MtpError("Cannot read alpha_moment_counts.");
    
    // alpha_index_basic_count
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_basic_count")
        MtpError("Cannot read alpha_index_basic_count.");
    ifs.ignore(2);
    ifs >> this->_alpha_index_basic_count;
    if (ifs.fail())
        MtpError("Cannot read alpha_index_basic_count.");
    
    // alpha_index_basic
    this->_alpha_index_basic = (int (*)[4])malloc(4*sizeof(int)*this->_alpha_index_basic_count);
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_basic")
        MtpError("Cannot read alpha_index_basic.");
    ifs.ignore(4);
    for (int ii=0; ii<this->_alpha_index_basic_count; ii++) {
        char tmpch;
        ifs.ignore(1000, '{');
        ifs >> this->_alpha_index_basic[ii][0] >> tmpch
            >> this->_alpha_index_basic[ii][1] >> tmpch
            >> this->_alpha_index_basic[ii][2] >> tmpch
            >> this->_alpha_index_basic[ii][3];
        if (ifs.fail())
            MtpError("Cannot read alpha_index_basic.");
    }
    ifs.ignore(1000, '\n');

    // alpha_index_times_count
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_times_count") 
        MtpError("Cannot read alpha_index_times_count.");
    ifs.ignore(2);
    ifs >> this->_alpha_index_times_count;
    if (ifs.fail())
        MtpError("Cannot read alpha_index_times_count.");

    // alpha_index_times
    if (this->_alpha_index_times_count != 0)
        this->_alpha_index_times = (int (*)[4])malloc(4*sizeof(int) * this->_alpha_index_times_count);
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_times")
        MtpError("Cannot read alpha_index_times.");
    ifs.ignore(4);
    if (this->_alpha_index_times_count != 0)
        for (int ii=0; ii<this->_alpha_index_times_count; ii++) {
            char tmpch;
            ifs.ignore(1000, '{');
            ifs >> this->_alpha_index_times[ii][0] >> tmpch
                >> this->_alpha_index_times[ii][1] >> tmpch
                >> this->_alpha_index_times[ii][2] >> tmpch
                >> this->_alpha_index_times[ii][3];
            if (ifs.fail())
                MtpError("Cannot read alpha_index_times.");
        }
    ifs.ignore(1000, '\n');

    // alpha_scalar_moments
    ifs >> tmpstr;
    if (tmpstr != "alpha_scalar_moments")
        MtpError("Cannot read alpha_scalar_moments.");
    ifs.ignore(2);
    ifs >> this->_alpha_scalar_moments;
    if (ifs.fail())
        MtpError("Cannot read alpha_scalar_moments.");
    
    // alpha_moment_mapping
    this->_alpha_moment_mapping = (int*)malloc(sizeof(int) * this->_alpha_scalar_moments);
    ifs >> tmpstr;
    if (tmpstr != "alpha_moment_mapping")
        MtpError("Cannot read alpha_moment_mapping.")
    ifs.ignore(4);
    for (int ii=0; ii<this->_alpha_scalar_moments; ii++) {
        char tmpch;
        ifs >> this->_alpha_moment_mapping[ii] >> tmpch;
        if (ifs.fail())
            MtpError("Cannot read alpha_moment_mapping.");
    }

    ifs.close();
}

MtpParam::MtpParam(const MtpParam& rhs)
{
    this->_alpha_moments_count = rhs._alpha_moments_count;
    this->_alpha_index_basic_count = rhs._alpha_index_basic_count;
    this->_alpha_index_basic = (int (*)[4])malloc(4*sizeof(int)*this->_alpha_index_basic_count);
    for (int ii=0; ii<this->_alpha_index_basic_count; ii++)
    {
        this->_alpha_index_basic[ii][0] = rhs._alpha_index_basic[ii][0];
        this->_alpha_index_basic[ii][1] = rhs._alpha_index_basic[ii][1];
        this->_alpha_index_basic[ii][2] = rhs._alpha_index_basic[ii][2];
        this->_alpha_index_basic[ii][3] = rhs._alpha_index_basic[ii][3];
    }
    this->_alpha_index_times_count = rhs._alpha_index_times_count;
    this->_alpha_index_times = (int (*)[4])malloc(4*sizeof(int)*this->_alpha_index_times_count);
    for (int ii=0; ii<this->_alpha_index_times_count; ii++)
    {
        this->_alpha_index_times[ii][0] = rhs._alpha_index_times[ii][0];
        this->_alpha_index_times[ii][1] = rhs._alpha_index_times[ii][1];
        this->_alpha_index_times[ii][2] = rhs._alpha_index_times[ii][2];
        this->_alpha_index_times[ii][3] = rhs._alpha_index_times[ii][3];
    }
    this->_alpha_scalar_moments = rhs._alpha_scalar_moments;
    this->_alpha_moment_mapping = (int*)malloc(sizeof(int) * this->_alpha_scalar_moments);
    for (int ii=0; ii<this->_alpha_scalar_moments; ii++)
        this->_alpha_moment_mapping[ii] = rhs._alpha_moment_mapping[ii];
}

MtpParam::MtpParam(MtpParam&& rhs)
{
    this->_alpha_moments_count = rhs._alpha_moments_count;
    this->_alpha_index_basic_count = rhs._alpha_index_basic_count;
    this->_alpha_index_basic = rhs._alpha_index_basic;
    this->_alpha_index_times_count = rhs._alpha_index_times_count;
    this->_alpha_index_times = rhs._alpha_index_times;
    this->_alpha_scalar_moments = rhs._alpha_scalar_moments;
    this->_alpha_moment_mapping = rhs._alpha_moment_mapping;


    rhs._alpha_moments_count = 0;
    rhs._alpha_index_basic_count = 0;
    rhs._alpha_index_basic = nullptr;
    rhs._alpha_index_times_count = 0;
    rhs._alpha_index_times = nullptr;
    rhs._alpha_scalar_moments = 0;
    rhs._alpha_moment_mapping = nullptr;
}

MtpParam& MtpParam::operator=(const MtpParam& rhs)
{
    this->_alpha_moments_count = 0;
    if (this->_alpha_index_basic_count != 0) {
        free(this->_alpha_index_basic);
        this->_alpha_index_basic_count = 0;
    }
    if (this->_alpha_index_times_count != 0) {
        free(this->_alpha_index_times);
        this->_alpha_index_times_count = 0;
    }
    if (this->_alpha_scalar_moments != 0) {
        free(this->_alpha_moment_mapping);
        this->_alpha_scalar_moments = 0;
    }

    this->_alpha_moments_count = rhs._alpha_moments_count;
    this->_alpha_index_basic_count = rhs._alpha_index_basic_count;
    if (this->_alpha_index_basic_count != 0) {
        this->_alpha_index_basic = (int (*)[4])malloc(4*sizeof(int)*this->_alpha_index_basic_count);
        for (int ii=0; ii<this->_alpha_index_basic_count; ii++) {
            this->_alpha_index_basic[ii][0] = rhs._alpha_index_basic[ii][0];
            this->_alpha_index_basic[ii][1] = rhs._alpha_index_basic[ii][1];
            this->_alpha_index_basic[ii][2] = rhs._alpha_index_basic[ii][2];
            this->_alpha_index_basic[ii][3] = rhs._alpha_index_basic[ii][3];
        }
    }
    this->_alpha_index_times_count = rhs._alpha_index_times_count;
    if (this->_alpha_index_times_count != 0) {
        this->_alpha_index_times = (int (*)[4])malloc(4*sizeof(int)*this->_alpha_index_times_count);
        for (int ii=0; ii<this->_alpha_index_times_count; ii++) {
            this->_alpha_index_times[ii][0] = rhs._alpha_index_times[ii][0];
            this->_alpha_index_times[ii][1] = rhs._alpha_index_times[ii][1];
            this->_alpha_index_times[ii][2] = rhs._alpha_index_times[ii][2];
            this->_alpha_index_times[ii][3] = rhs._alpha_index_times[ii][3];
        }
    }
    this->_alpha_scalar_moments = rhs._alpha_scalar_moments;
    if (this->_alpha_scalar_moments != 0) {
        this->_alpha_moment_mapping = (int*)malloc(sizeof(int) * this->_alpha_scalar_moments);
        for (int ii=0; ii<this->_alpha_scalar_moments; ii++)
            this->_alpha_moment_mapping[ii] = rhs._alpha_moment_mapping[ii];
    }

    return *this;
}

MtpParam& MtpParam::operator=(MtpParam&& rhs)
{
    if (this != &rhs) {
        this->_alpha_moments_count = 0;
        if (this->_alpha_index_basic_count != 0) {
            free(this->_alpha_index_basic);
            this->_alpha_index_basic_count = 0;
        }
        if (this->_alpha_index_times_count != 0) {
            free(this->_alpha_index_times);
            this->_alpha_index_times_count = 0;
        }
        if (this->_alpha_scalar_moments != 0) {
            free(this->_alpha_moment_mapping);
            this->_alpha_scalar_moments = 0;
        }

        this->_alpha_moments_count = rhs._alpha_moments_count;
        this->_alpha_index_basic_count = rhs._alpha_index_basic_count;
        this->_alpha_index_basic = rhs._alpha_index_basic;
        this->_alpha_index_times_count = rhs._alpha_index_times_count;
        this->_alpha_index_times = rhs._alpha_index_basic;
        this->_alpha_scalar_moments = rhs._alpha_scalar_moments;
        this->_alpha_moment_mapping = rhs._alpha_moment_mapping;

        rhs._alpha_moments_count = 0;
        rhs._alpha_index_basic_count = 0;
        rhs._alpha_index_basic = nullptr;
        rhs._alpha_index_times_count = 0;
        rhs._alpha_index_times = nullptr;
        rhs._alpha_scalar_moments = 0;
        rhs._alpha_moment_mapping = nullptr;
    }

}

MtpParam::~MtpParam()
{
    free(this->_alpha_index_basic);
    free(this->_alpha_index_times);
    free(this->_alpha_moment_mapping);
}

void MtpParam::show() const
{
    printf("1. alpha_moments_count = %10d\n", this->_alpha_moments_count);
    printf("2. alpha_index_basic_count = %10d\n", this->_alpha_index_basic_count);
    printf("3. alpha_index_basic:\n");
    for (int ii=0; ii<this->_alpha_index_basic_count; ii++)
        printf("\t[%5d, %5d, %5d, %5d]\n", this->_alpha_index_basic[ii][0], this->_alpha_index_basic[ii][1], this->_alpha_index_basic[ii][2], this->_alpha_index_basic[ii][3]);
    printf("4. alpha_index_times_count = %10d\n", this->_alpha_index_times_count);
    printf("5. alpha_index_times:\n");
    for (int ii=0; ii<this->_alpha_index_times_count; ii++)
        printf("\t[%5d, %5d, %5d, %5d]\n", this->_alpha_index_times[ii][0], this->_alpha_index_times[ii][1], this->_alpha_index_times[ii][2], this->_alpha_index_times[ii][3]);
    printf("6. alpha_scalar_moments = %10d\n", this->_alpha_scalar_moments);
    printf("7. alpha_moment_mapping:\n\t");
    for (int ii=0; ii<this->_alpha_scalar_moments; ii++)
        printf("%5d, ", this->_alpha_moment_mapping[ii]);
    printf("\n");
}

const int MtpParam::alpha_moments_count() const
{
    return this->_alpha_moments_count;
}

const int MtpParam::alpha_index_basic_count() const
{
    return this->_alpha_index_basic_count;
}

const int (*MtpParam::alpha_index_basic() const)[4]
{
    return this->_alpha_index_basic;
}

const int MtpParam::alpha_index_times_count() const
{
    return this->_alpha_index_times_count;
}

const int (*MtpParam::alpha_index_times() const)[4]
{
    return this->_alpha_index_times;
}

const int MtpParam::alpha_scalar_moments() const 
{
    return this->_alpha_scalar_moments;
}

const int *MtpParam::alpha_moment_mapping() const
{
    return this->_alpha_moment_mapping;
}

};  // namespace : mtpr
};  // namespace : matersdk