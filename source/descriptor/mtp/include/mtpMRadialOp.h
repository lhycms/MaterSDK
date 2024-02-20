#ifndef MATERSDK_MTPM_RADIAL_OP
#define MATERSDK_MTPM_RADIAL_OP
#include <torch/torch.h>
#include "./mtpMRadial.h"

namespace matersdk {
namespace mtp {

class MtpQFunction : public torch::autograd::Function<MtpQFunction>
{
public:
    static torch::autograd::variable_list forward(
        torch::autograd::AutogradContext* ctx,
        int64_t size,
        at::Tensor rs_tensor);  // rs_tensor = [rcut, rcut_smooth, distance_ij]
    
    static torch::autograd::variable_list backward(
        torch::autograd::AutogradContext* ctx,
        torch::autograd::variable_list grad_outputs);
};  // class MtpQFunction


torch::autograd::variable_list MtpQOp(
    int64_t size,
    at::Tensor rs_tensor);

};  // namespace : mtp
};  // namespace : matersdk

#endif