import torch
import torch.nn as nn
from typing import List


class DpEmbeddingNet(nn.Module):
    def __init__(
                self,
                input_dims:List[int],
                embedding_sizes:List[int],
                M2:int,
                acticvation_fn=nn.Tanh):
        '''
        Description
        ----------- 
            1. DeepPot-SE: Embedding net implementation
            2. The shape of descriptor of DeepPot-SE is:
                - (num_frames, max_num_nbrs, M1, M2)
                - `M1 = embedding_sizes[-1]`
                - `M2 << M1`
        
        Parameters
        ----------
            1. input_dim: List[int]
                - If using `s(r_ji)`, `input_size = [0]`
                - If using `x_ji, y_ji, z_ji`, `input_size = [1, 2, 3]`
                - If using `s(r_ji), x_ji, y_ji, z_ji`, `input_size = [0, 1, 2, 3]`
            2. embedding_sizes: List[int]
                - Hidden layer of embedding net
            3. M2: int
                - 
        '''
        super(DpEmbeddingNet, self).__init__()
        ### Step 1. 初始化自身属性
        self.input_dims = input_dims
        self.sizes = [len(input_dims)] + embedding_sizes
        self.M2 = M2
        self.activation_fn = acticvation_fn()
        self.weights = nn.ParameterDict()
        self.bias = nn.ParameterDict()
        
        ### Step 2. 用 `nn.init.normal_(torch.Tensor, mean:float, std:float)` 初始化 wij, bi
        for tmp_i in range(1, len(self.sizes)):
            ### Step 2.1. Weights
            tmp_weights = torch.Tensor(self.sizes[tmp_i-1], self.sizes[tmp_i])
            nn.init.normal_(tmp_weights, mean=0, std=0.01)
            self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_weights)
            
            ### Step 2.2. Bias
            tmp_bias = torch.Tensor(self.sizes[tmp_i])
            nn.init.normal_(tmp_bias, mean=0, std=0.01)
            self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_bias)
    
    
    def forward(self, tilde_r:torch.Tensor):
        '''
        Description
        -----------
            1. You will get G from `self.forward()`
            2. You can use `G` and `R` to calculate `Descriptor - D`, then use it in fitting net.
                - D^i = G^i^T \cdot R \cdot R^T \cdot G^i
        '''
        num_frames, num_centers, max_num_nbrs, num_r = tilde_r.size()
        x = tilde_r[:, :, :, self.input_dims].view(num_frames, num_centers, max_num_nbrs, len(self.input_dims))
        print("Size = ", x.size())
        
        for tmp_i in range(1, len(self.sizes)):
            x = torch.matmul(x, self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)])
            x = x + self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)].view(1, -1)
            x = self.activation_fn(x)
        G = x
        
        #D = torch.matmul(G.transpose(-1, -2), tilde_r)
        #D = torch.matmul(D, tilde_r.transpose(-1, -2))
        #D = torch.matmul(D, G[:, :, :, :self.M2])
        
        return G
    
    
    def check(self):
        for tmp_i in range(1, len(self.sizes)):
            print('*' * 40)
            print("w_{0}_{1} = ".format(tmp_i-1, tmp_i), self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)].size())
            print("b_{0}_{1} = ".format(tmp_i-1, tmp_i), self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)].size())
            print('*' * 40)