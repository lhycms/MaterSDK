#include <torch/torch.h>
#include "../include/mtpMRadial.h"
#include "../include/mtpMRadialOp.h"


namespace matersdk {
namespace mtp {

torch::autograd::variable_list MtpQFunction::forward(
    torch::autograd::AutogradContext* ctx,
    int64_t size,
    at::Tensor rs_tensor)
{
    c10::TensorOptions options = c10::TensorOptions()
        .dtype(rs_tensor.scalar_type())
        .device(rs_tensor.device());
    at::Tensor result_tensor = at::zeros({size}, options);
    at::Tensor deriv2r_tensor = at::zeros({size}, options);
    if (options.dtype() == torch::kFloat32) {
        MtpQ<float> mtp_q(
            size, 
            rs_tensor[0].item<float>(), 
            rs_tensor[1].item<float>());
        mtp_q.build(rs_tensor[2].item<float>());
        
        float* result = result_tensor.data_ptr<float>();
        float* deriv2r = deriv2r_tensor.data_ptr<float>();
        for (int ii=0; ii<(int)size; ii++) {
            result[ii] = mtp_q.get_result()[ii];
            deriv2r[ii] = mtp_q.get_deriv2r()[ii];
        }
    } else {
        MtpQ<double> mtp_q(
            size,
            rs_tensor[0].item<double>(),
            rs_tensor[1].item<double>());
        mtp_q.build(rs_tensor[2].item<double>());
        double* result = result_tensor.data_ptr<double>();
        double* deriv2r = deriv2r_tensor.data_ptr<double>();
        for (int ii=0; ii<(int)size; ii++) {
            result[ii] = mtp_q.get_result()[ii];
            deriv2r[ii] = mtp_q.get_deriv2r()[ii];
        }
    }

    ctx->save_for_backward({deriv2r_tensor});
    return {result_tensor};
}


torch::autograd::variable_list MtpQFunction::backward(
    torch::autograd::AutogradContext* ctx,
    torch::autograd::variable_list grad_outputs)
{
    at::Tensor deriv2r_tensor = ctx->get_saved_variables()[0];
    at::Tensor grad_output_tensor = grad_outputs[0];
    if (!grad_output_tensor.is_contiguous()) {
        grad_output_tensor = grad_output_tensor.contiguous();
    }
    assert(deriv2r_tensor.scalar_type() == grad_output_tensor.scalar_type());

    c10::TensorOptions options = c10::TensorOptions()
        .dtype(deriv2r_tensor.scalar_type())
        .device(deriv2r_tensor.device());
    at::Tensor deriv_tensor = at::zeros({1}, options);
    if (options.dtype() == torch::kFloat32) {
        float* deriv2r = deriv2r_tensor.data_ptr<float>();
        float* grad_output = grad_output_tensor.data_ptr<float>();
        float* deriv = deriv_tensor.data_ptr<float>();
        for (int ii=0; ii<(int)deriv_tensor.sizes()[0]; ii++) 
            deriv[0] = grad_output[ii] * deriv2r[ii];
    } else {
        double* deriv2r = deriv2r_tensor.data_ptr<double>();
        double* grad_output = grad_output_tensor.data_ptr<double>();
        double* deriv = deriv_tensor.data_ptr<double>();
        for (int ii=0; ii<(int)deriv_tensor.sizes()[0]; ii++)
            deriv[0] = grad_output[ii] * deriv2r[ii];
    }
    return {deriv_tensor};
}


torch::autograd::variable_list MtpQOp(
    int64_t size,
    at::Tensor rs_tensor)
{
    return MtpQFunction::apply(size, rs_tensor);
}

};  // class : mtp
};  // class : matersdk