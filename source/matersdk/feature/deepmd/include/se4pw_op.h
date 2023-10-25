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
        at::Tensor& types,  // x, types 一一对应
        int ntypes,
        at::Tensor& num_neigh_atoms_lst,
        double rcut,             // Python 中默认是双精度浮点数
        double rcut_smooth);     // Python 中默认是双精度浮点数
};

};  // namespace : deepPotSE
};  // namespace : matersdk


#endif