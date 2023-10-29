#include <torch/torch.h>
#include <cstdlib>
#include <stdio.h>
#include "../include/se4pw_op.h"


namespace matersdk {
namespace deepPotSE {


/**
 * @brief 
 * 
 * @note It doesn't matter whether input tensor is flatten or not.
 *  1. int: torch::kInt32
 *  2. float: torch::kFloat32, torch::kFloat64
 * 
 * @param inum 
 * @param ilist 
 * @param numneigh 
 * @param firstneigh 
 * @param x 
 * @param types 
 * @param ntypes 
 * @param num_neigh_atoms_lst 
 * @param rcut 
 * @param rcut_smooth 
 * @return torch::autograd::variable_list 
 */
torch::autograd::variable_list Se4pwOp::forward(
        int inum,
        at::Tensor& ilist,
        at::Tensor& numneigh,
        at::Tensor& firstneigh,
        at::Tensor& x,
        at::Tensor& types,
        int ntypes,
        at::Tensor& num_neigh_atoms_lst,
        double rcut,
        double rcut_smooth)
{
    c10::Device device = x.device();
    c10::ScalarType dtype = x.scalar_type();
    c10::TensorOptions tensor_options = c10::TensorOptions().device(device).dtype(dtype);
    c10::ScalarType int_dtype = ilist.scalar_type();
    
    ilist.to(int_dtype);
    numneigh.to(int_dtype);
    firstneigh.to(int_dtype);
    types.to(int_dtype);
    num_neigh_atoms_lst.to(int_dtype);
    
    int tot_num_neigh_atoms = num_neigh_atoms_lst.sum().item<int>();
    at::Tensor tilde_r = at::zeros({inum, tot_num_neigh_atoms, 4}, tensor_options);
    at::Tensor tilde_r_deriv = at::zeros({inum, tot_num_neigh_atoms, 4, 3}, tensor_options);
    at::Tensor relative_coords = at::zeros({inum, tot_num_neigh_atoms, 3}, tensor_options);

    if (dtype == torch::kFloat32) {
        float* tilde_r_ptr = tilde_r.data_ptr<float>();
        float* tilde_r_deriv_ptr = tilde_r_deriv.data_ptr<float>();
        float* relative_coords_ptr = relative_coords.data_ptr<float>();

        Se4pw<float>::generate(
                tilde_r_ptr,
                tilde_r_deriv_ptr,
                relative_coords_ptr,
                inum,
                ilist.data_ptr<int>(),
                numneigh.data_ptr<int>(),
                firstneigh.data_ptr<int>(),
                x.data_ptr<float>(),
                types.data_ptr<int>(),
                ntypes,
                num_neigh_atoms_lst.data_ptr<int>(),
                (float)rcut,
                (float)rcut_smooth);
    } else {
        double* tilde_r_ptr = tilde_r.data_ptr<double>();
        double* tilde_r_deriv_ptr = tilde_r_deriv.data_ptr<double>();
        double* relative_coords_ptr = relative_coords.data_ptr<double>();
        Se4pw<double>::generate(
                tilde_r_ptr,
                tilde_r_deriv_ptr,
                relative_coords_ptr,
                inum,
                ilist.data_ptr<int>(),
                numneigh.data_ptr<int>(),
                firstneigh.data_ptr<int>(),
                x.data_ptr<double>(),
                types.data_ptr<int>(),
                ntypes,
                num_neigh_atoms_lst.data_ptr<int>(),
                rcut,
                rcut_smooth);
    }


    return {tilde_r, tilde_r_deriv, relative_coords};
}


};  // namespace : deepPotSE 
};  // namespace : matersdk