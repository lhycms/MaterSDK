#ifndef MATERSDK_MTPR_MTP_PARAM_H
#define MATERSDK_MTPR_MTP_PARAM_H
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

    MtpParam(const MtpParam& rhs);

    MtpParam(MtpParam&& rhs);

    MtpParam& operator=(const MtpParam& rhs);

    MtpParam& operator=(MtpParam&& rhs);

    ~MtpParam();

    void show() const;

    const int alpha_moments_count() const;

    const int alpha_index_basic_count() const;

    const int (*alpha_index_basic() const)[4];

    const int alpha_index_times_count() const;

    const int (*alpha_index_times() const)[4];

    const int alpha_scalar_moments() const;

    const int *alpha_moment_mapping() const;

    const int nmus() const;

private:
    int _alpha_moments_count = 0;
    int _alpha_index_basic_count = 0;
    int (*_alpha_index_basic)[4] = nullptr;
    int _alpha_index_times_count = 0;
    int (*_alpha_index_times)[4] = nullptr;
    int _alpha_scalar_moments = 0;
    int *_alpha_moment_mapping = nullptr;
    //int _alpha_count = 0;   // Basis function count. _alpha_count = _alpha_scalar_moments + 1
};  // class MtpParam


class AlphaIndexTimes {
public:
    AlphaIndexTimes();

    AlphaIndexTimes(
        int alpha_index_times_count,
        int (*alpha_index_times)[4]);
    
    AlphaIndexTimes(const AlphaIndexTimes &rhs);

    AlphaIndexTimes(AlphaIndexTimes &&rhs);

    AlphaIndexTimes& operator=(const AlphaIndexTimes &rhs);

    AlphaIndexTimes& operator=(AlphaIndexTimes &&rhs);

    ~AlphaIndexTimes();

    void show() const;

    const int alpha_index_times_count() const;

    const int (*alpha_index_times() const)[4];

private:
    int _alpha_index_times_count = 0;
    int (*_alpha_index_times)[4] = nullptr;
};  // class : AlphaIndexTimes

std::set<int> find_mus4nonbasic_mom(
    const int mom_idx,
    const int alpha_moments_count,
    const int alpha_index_basic_count,
    const int (*alpha_index_basic)[4],
    const int alpha_index_times_count,
    const int (*alpha_index_times)[4],
    const int alpha_scalar_moments,
    const int *alpha_moment_mapping);


};  // namespace : mtpr
};  // namespace : matersdk

#endif