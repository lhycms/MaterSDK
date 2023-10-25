#ifndef MATERSDK_SE4PW_OP_H
#define MATERSDK_SE4PW_OP_H

#include "./se4pw.h"
#include <torch/torch.h>
#include <cstdlib>

namespace matersdk {
namespace deepPotSE {
class Se4pwOp {
public:
    static torch::autograd::variable_list forward(
        int inum,
        at::Tensor& ilist,
        at::Tensor& numneigh,
        at::Tensor& firstneigh,
        at::Tensor& x,
        at::Tensor& types,
        int ntypes,
        at::Tensor& num_neigh_atoms_lst,
        double rcut,            // Python 中默认是双精度浮点数
        double rcut_smooth)     // Python 中默认是双精度浮点数
    {
        int tot_num_neigh_atoms = 0;
        for (int ii=0; ii<ntypes; ii++)
            tot_num_neigh_atoms += num_neigh_atoms_lst[ii]; // 特征的维度
        int num_cols_numneigh = firstneigh.sizes()[1];      // firstneigh 的 sizes()[1]，不是实际最大近邻数，因为firstneigh后面可能用 -1 填充
        // Step 1. determine dtype and device, init `tilde_r`, `tilde_r_deriv`, `relative_coords`
        at::ScalarType dtype = x.scalar_type();
        at::ScalarType dtype_int = num_neigh_atoms_lst.scalar_type();
        at::Device device = x.device();
        c10::TensorOptions tensor_options = c10::TensorOptions().dtype(dtype).device(device);
        at::Tensor tilde_r = at::Tensor({inum, tot_num_neigh_atoms, 4}, tensor_options);
        at::Tensor tilde_r_deriv = at::Tensor({inum, tot_num_neigh_atoms, 4, 3}, tensor_options);
        at::Tensor relative_coords = at::Tensor({inum, tot_num_neigh_atoms, 3}, tensor_options);
        tilde_r.reuiqres_grad_(true);
        tilde_r_deriv.requires_grad_(true);
        relative_coords.requires_grad_(true);

        // Step 2. Prepare input for pure C++ code (matersdk::deepPotSE::Se4pw)
        // Step 2.1. Get pointer to data of at::Tensor
        if (dtype == torch::kFloat32) {
            float* tilde_r_ptr = tilde_r.data_ptr<float>();
            float* tilde_r_deriv_ptr = tilde_r_deriv.data_ptr<float>();
            float* relative_coords_ptr = relative_coords.data_ptr<float>();
            float* x_ptr = x.data_ptr<float>();
        } else (dtype == torch::kFloat64) {
            double* tilde_r_ptr = tilde_r.data_ptr<double>();
            double* tilde_r_deriv_ptr = tilde_r_deriv.data_ptr<double>();
            double* relative_coords_ptr = relative_coords.data_ptr<double>();
            double* x_ptr = x.data_ptr<double>();
        }

        if (dtype_int == torch::kInt16) {
            int16_t* ilist_ptr = ilist.data_ptr<int16_t>();
            int16_t* numneigh_ptr = numneigh.data_ptr<int16_t>();
            int16_t* firstneigh_ptr = firstneigh.data_ptr<int16_t>();
            int16_t* types_ptr = types.data_ptr<int16_t>();
            int16_t* num_neigh_atoms_lst_ptr = num_neigh_atoms_lst.data_ptr<int16_t>();
            
            int16_t** firstneigh_2dptr = (int16_t**)malloc(sizeof(int16_t*) * inum);
            for (int ii=0; ii<inum; ii++)
                firstneigh_2dptr[ii] = (int16_t*)malloc(sizeof(int16_t) * max_numneigh);
            for (int ii=0; ii<inum; ii++)
                for (int jj=0; jj<numneigh[ii]; jj++)
                    firstneigh_2dptr[ii][jj] = firstneigh_ptr[ii*num_cols_numneigh+jj];

        } else if (dtype_int == torch::kInt32) {
            int32_t* ilist_ptr = ilist.data_ptr<int32_t>();
            int32_t* numneigh_ptr = numneigh.data_ptr<int32_t>();
            int32_t* firstneigh_ptr = firstneigh.data_ptr<int32_t>();
            int32_t* types_ptr = types.data_ptr<int32_t>();
            int32_t* num_neigh_atoms_lst_ptr = num_neigh_atoms_lst.data_ptr<int32_t>();
        } else (dtype_int == torch::kInt64) {
            int64_t* ilist_ptr = ilist.data_ptr<int64_t>();
            int64_t* numneigh_ptr = numneigh.data_ptr<int64_t>();
            int64_t* firstneigh_ptr = firstneigh.data_ptr<int64_t>();
            int64_t* types_ptr = types.data_ptr<int64_t>();
            int64_t* num_neigh_atoms_lst_ptr = num_neigh_atoms_lst.data_ptr<int64_t>();
        }

        // Step 2.2. Convert 1d pointer to 2d pointer (firstneigh, x)
        

        // Step 3. 


        // Step . Free memory
    }
};

};  // namespace : deepPotSE
};  // namespace : matersdk