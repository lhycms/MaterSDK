#ifndef MATERSDK_MTP_PARAM_H
#define MATERSDK_MTP_PARAM_H
#include <string>
#include <vector>

namespace matersdk {
namespace mtpr {

class MtpException {
public:
    std::string message;

    MtpException(const std::string& message) : message(message) {}

    virtual const char* what() const {
        return this->message.c_str();
    } 
};  // class : MtpException

#define MtpError(str) throw MtpException((std::string)"ERROR: " + str + \
    " \n    Thrown by function " + __FUNCTION__ + \
    ",\n    source file " + __FILE__ + \
    ",\n    line " + std::to_string(static_cast<long long>(__LINE__)) + "\n");

class MtpParam {
public:
    MtpParam();

    MtpParam(const std::string& filename);

    void _load(const std::string& filename);

    ~MtpParam();

private:
    std::string pot_desc;
    double scaling = 1.0;
    int species_count = 0;
    std::string rbasis_type;
    std::string mbasis_type;
    int radial_func_count;
    int _alpha_moments_count = 0;
    int _alpha_index_basic_count = 0;
    int (*_alpha_index_basic)[4];
    int _alpha_index_times_count = 0;
    int (*_alpha_index_times)[4];
    int _alpha_scalar_moments = 0;
    int* _alpha_moment_mapping;
    int _alpha_count = 0;   // Basis function count. _alpha_count = _alpha_scalar_moments + 1
};  // class MtpParam

};  // namespace : mtpr
};  // namespace : matersdk

#endif