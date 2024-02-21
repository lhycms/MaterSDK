#ifndef MATERSDK_MTPF_MODULE_H
#define MATERSDK_MTPF_MODULE_H
#include <torch/torch.h>

namespace matersdk {
namespace mtp {
class MtpFModuleImpl : public torch::nn::Module
{
public:
    MtpFModuleImpl(
        int64_t nmus,
        int64_t ntypes,
        int64_t size,
        at::Tensor rcuts_tensor);
    
    at::Tensor forward(
        int64_t mu,
        int64_t iidx,
        at::Tensor firstneigh_tensor,
        at::Tensor types,
        at::Tensor rcs_tensor);

private:
    int64_t nmus;
    int64_t ntypes;
    int64_t size;
    at::Tensor rcuts_tensor;
    torch::nn::ParameterList cheby_coeff_list;
};  // class : MtpFModuleImpl

TORCH_MODULE(MtpFModule);

};  // namespace : matersdk
};  // namespace : mtp

#endif