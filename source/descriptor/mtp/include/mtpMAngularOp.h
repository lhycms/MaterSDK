#ifndef MATERSDK_MTPM_ANGULAR_OP_H
#define MATERSDK_MTPM_ANGULAR_OP_H
#include <torch/torch.h>

namespace matersdk {
namespace mtp {

at::Tensor MtpMAngularOp(
    const at::Tensor& relative_coord_tensor,
    int nu);

}; // namespace : mtp
}; // namespace : matersdk

#endif
