#include <torch/torch.h>
#include "../include/mtpMRadialOp.h"
#include "../include/mtpFModule.h"


namespace matersdk {
namespace mtp {

MtpFModuleImpl::MtpFModuleImpl(
    int64_t nmus,
    int64_t ntypes,
    int64_t size,
    at::Tensor rcuts_tensor)
{
    this->size = size;
    this->nmus = nmus;
    this->ntypes = ntypes;
    this->rcuts_tensor = rcuts_tensor;
    int64_t ncoeffs = nmus * ntypes * ntypes;
    for (int64_t ii=0; ii<ncoeffs; ii++)
        this->cheby_coeff_list->append(at::randn({size}));
    this->register_module("cheby_coeff_list", this->cheby_coeff_list);
}

/**
 * @brief Calculate f_{\mu}(|r_{ij}, z_i, z_j|)
 * @param iidx The index of center atoms
 * @param ifirstneigh_tensor The indices of neigh atoms.
 *          .shape = [nneighs,]
 * @param types The atomic types
 * @param ircs_tensor The relative coordinates of neigh atoms
 *          .shape = [nneighs, 3]
 * @return at::Tensor 
 *          .shape = [nneighs,]
 */
at::Tensor MtpFModuleImpl::forward(
    int64_t mu,
    int64_t iidx,
    at::Tensor ifirstneigh_tensor,
    at::Tensor types,
    at::Tensor ircs_tensor)
{
    at::Tensor mtp_f_tensor = at::zeros({ifirstneigh_tensor.sizes()[0]})
        .requires_grad_(false);
    int64_t z_i = types[iidx].item<int64_t>();
    int64_t* firstneigh = ifirstneigh_tensor.data_ptr<int64_t>();    
    at::Tensor mtp_q_tensor = MtpQOp(   // .shape = [nneigh, size]
        this->size,
        this->rcuts_tensor,
        ircs_tensor)[0];
std::cout << mtp_q_tensor << std::endl;
    for (int64_t ii=0; ii<ifirstneigh_tensor.sizes()[0]; ii++) {
        int64_t z_j = types[firstneigh[ii]].item<int64_t>();
        mtp_f_tensor[ii] = at::dot(
            this->cheby_coeff_list[mu*this->ntypes*this->ntypes + z_i*this->ntypes + z_j],
            mtp_q_tensor[ii]
        );
    }
    return mtp_f_tensor;
}

};  // namespace : mtp
};  // namespace : matersdk