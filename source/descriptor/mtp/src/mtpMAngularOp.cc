#include <cassert>
#include <cmath>
#include "../include/mtpMAngularOp.h"

namespace matersdk {
namespace mtp {


/**
 * @brief Given i,j, and \nu. Calculate $r_{ij} \otimes_{nu} r_{ij}$
 * 
 * @param relative_coord_tensor sizes() = (3)
 * @param nu $\nu$ in r_{ij} \otimes_{nu} r_{ij}
 * @return at::Tensor 
 */
at::Tensor MtpMAngularOp(
    const at::Tensor& relative_coord_tensor,
    int nu)
{
    c10::TensorOptions options = c10::TensorOptions()
        .device(relative_coord_tensor.device())
        .dtype(relative_coord_tensor.scalar_type());
    if (nu == 0) {
        return at::ones({1}, options);
    } else if (nu == 1) {
        return relative_coord_tensor;
    } else {
        at::Tensor mtp_angular = relative_coord_tensor.clone();
        for (int ii=0; ii<nu-1; ii++)
            mtp_angular = torch::tensordot(
                mtp_angular, 
                relative_coord_tensor,
                {}, {});
        return mtp_angular;
    }
}

};  // namespace : mtp
};  // namespace : matersdk