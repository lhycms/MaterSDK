#ifndef MATERSDK_ANGULAR_OP_H
#define MATERSDK_ANGULAR_OP_H
#include <torch/torch.h>

namespace matersdk {
namespace mtp {

at::Tensor MtpMAngularOp(
    at::Tensor relative_coord_tensor,
    int nu);

}; // namespace : mtp
}; // namespace : matersdk

#endif