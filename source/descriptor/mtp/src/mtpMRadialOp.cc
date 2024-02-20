#include <torch/torch.h>
#include <vector>
#include <iterator>
#include "../include/mtpMRadial.h"
#include "../include/mtpMRadialOp.h"


namespace matersdk {
namespace mtp {

torch::autograd::variable_list MtpQFunction::forward(
    torch::autograd::AutogradContext* ctx,
    int64_t size,
    at::Tensor rcuts_tensor,
    at::Tensor distances_tensor)
{
    c10::TensorOptions options = c10::TensorOptions()
        .dtype(distances_tensor.scalar_type())
        .device(distances_tensor.device());
    at::IntArrayRef bcn_shape = distances_tensor.sizes();   // [nneighs,]
    std::vector<int64_t> mtp_q_shape;
    mtp_q_shape.clear();
    std::copy(bcn_shape.begin(), bcn_shape.end(), std::back_inserter<std::vector<int64_t>>(mtp_q_shape));
    mtp_q_shape.push_back(size);    // [nneighs, size]
    at::Tensor result_tensor = at::zeros(mtp_q_shape, options);
    at::Tensor deriv2r_tensor = at::zeros(mtp_q_shape, options);
    if (options.dtype() == torch::kFloat32) {
        float* distances = distances_tensor.data_ptr<float>();
        MtpQ<float> mtp_q(
            size, 
            rcuts_tensor[0].item<float>(), 
            rcuts_tensor[1].item<float>());
        float* result = result_tensor.data_ptr<float>();
        float* deriv2r = deriv2r_tensor.data_ptr<float>();
        for (int ii=0; ii<(int)bcn_shape[0]; ii++) {
            mtp_q.build(distances[ii]);
            for (int jj=0; jj<size; jj++) {
                result[ii*size + jj] = mtp_q.get_result()[jj];
                deriv2r[ii*size + jj] = mtp_q.get_deriv2r()[jj];
            }
        }
        
    } else {
        double* distances = distances_tensor.data_ptr<double>();
        MtpQ<double> mtp_q(
            size,
            rcuts_tensor[0].item<double>(),
            rcuts_tensor[1].item<double>());
        double* result = result_tensor.data_ptr<double>();
        double* deriv2r = deriv2r_tensor.data_ptr<double>();
        for (int ii=0; ii<(int)bcn_shape[0]; ii++) {
            mtp_q.build(distances[ii]);
            for (int jj=0; jj<size; jj++) {
                result[ii*size + jj] = mtp_q.get_result()[jj];
                deriv2r[ii*size + jj] = mtp_q.get_deriv2r()[jj];
            }
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
    at::Tensor deriv_tensor = at::zeros({deriv2r_tensor.sizes()[0]}, options);
    int64_t num_neighs = deriv2r_tensor.sizes()[0]; // nneighs
    int64_t size = deriv2r_tensor.sizes()[1];       // size
    if (options.dtype() == torch::kFloat32) {
        float* deriv2r = deriv2r_tensor.data_ptr<float>();
        float* grad_output = grad_output_tensor.data_ptr<float>();
        float* deriv = deriv_tensor.data_ptr<float>();
        for (int ii=0; ii<(int)num_neighs; ii++) 
            for (int jj=0; jj<(int)size; jj++)
                deriv[ii] += grad_output[ii*size + jj] * deriv2r[ii*size + jj];
    } else {
        double* deriv2r = deriv2r_tensor.data_ptr<double>();
        double* grad_output = grad_output_tensor.data_ptr<double>();
        double* deriv = deriv_tensor.data_ptr<double>();
std::cout << deriv2r_tensor << std::endl;
        for (int ii=0; ii<(int)num_neighs; ii++) 
            for (int jj=0; jj<(int)size; jj++)
                deriv[ii] += grad_output[ii*size + jj] * deriv2r[ii*size + jj];
    }
    return {at::Tensor(), at::Tensor(), deriv_tensor};
}


torch::autograd::variable_list MtpQOp(
    int64_t size,
    at::Tensor rcuts_tensor,
    at::Tensor distances_tensor)
{
    return MtpQFunction::apply(size, rcuts_tensor, distances_tensor);
}

};  // class : mtp
};  // class : matersdk
