import torch
import torch.nn as nn
from typing import List


class DpEmbeddingNet(nn.Module):
    def __init__(
                self,
                input_size:int,
                embedding_sizes:List[int],
                M2:int,
                activation_fn=nn.Tanh(),
                resnet_dt:bool=False):
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
            1. input_size: int
                - If using `s(r_ji)`, `input_size = 1`
                - If using `x_ji, y_ji, z_ji`, `input_size = 3`
                - If using `s(r_ji), x_ji, y_ji, z_ji`, `input_size = 4`
            2. embedding_sizes: List[int]
                - Hidden layer of embedding net
            3. M2: int
                - 
            4. activation_fn: nn.Tanh()
                - 激活函数
            5. resnet_dt: bool
                - 是否采用 `residual connection`: y=f(x)+x
        '''
        super(DpEmbeddingNet, self).__init__()
        
        ### Step 1. 初始化 EmbeddingNet 的属性
        self.sizes = [input_size] + embedding_sizes
        self.M2 = M2
        self.activation_fn = activation_fn
        self.resnet_dt = resnet_dt
        
        ### Step 2. 初始化 EmbeddingNet 的网络参数
        self.weights = nn.ParameterDict()
        self.bias = nn.ParameterDict()
        self.resnet_dts = nn.ParameterDict()
        
        for tmp_i in range(1, len(self.sizes)):
            ### Step 2.1. weight
            tmp_weights = torch.Tensor(self.sizes[tmp_i-1], self.sizes[tmp_i])
            nn.init.normal_(tmp_weights, mean=0, std=0.01)
            self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_weights, requires_grad=True)
            
            ### Step 2.2. bias
            tmp_bias = torch.Tensor(self.sizes[tmp_i])
            nn.init.normal_(tmp_bias, mean=0, std=0.01)
            self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_bias, requires_grad=True)
            
            ### Step 2.3. resnet_dts
            tmp_resnet_dts = torch.Tensor(self.sizes[tmp_i])
            nn.init.normal_(tmp_resnet_dts, mean=0, std=0.01)
            self.resnet_dts["r_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_resnet_dts, requires_grad=True)
            
    
    def forward(self, xx:torch.Tensor):
        '''
        Description
        -----------
            1. You will get G from `self.forward()`
            2. You can use `G` and `R` to calculate `Descriptor - D`, then use it in fitting net.
                - D^i = G^i^T \cdot R \cdot R^T \cdot G^i
        
        Parameters
        ----------
            1. xx(tilde_r_part): torch.Tensor
                - If using `s(r_ji)`:                   tilde_r_part.shape = (num_frames, num_centers, max_num_nbrs, 1)
                - If using `x_ji, y_ji, z_ji`:          tilde_r_part.shape = (num_frames, num_centers, max_num_nbrs, 3)    
                - If using `s(r_ji), x_ji, y_ji, z_ji`: tilde_r_part.shape = (num_frames, num_centers, max_num_nbrs, 4)
        '''
        assert (xx.size()[-1] == self.sizes[0])
        
        for tmp_i in range(1, len(self.sizes)):
            ### Step 1. Linear operation
            hidden = torch.matmul(xx, self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)])
            hidden = hidden + self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)].view(1, -1)
            
            ### Step 2. Activation function
            hidden = self.activation_fn(hidden)
            
            ### Step 3. Resnet block
            if self.sizes[tmp_i] == self.sizes[tmp_i-1]:
                if self.resnet_dt:
                    xx = hidden * self.resnet_dts['r_{0}_{1}'.format(tmp_i-1, tmp_i)].view(1, -1) + xx
                else:
                    xx = hidden + xx
            elif self.sizes[tmp_i] == 2 * self.sizes[tmp_i-1]:
                if self.resnet_dt:
                    xx = hidden * self.resnet_dts['r_{0}_{1}'.format(tmp_i-1, tmp_i)].view(1, -1) + torch.cat((xx, xx), dim=-1)
                else:
                    xx = hidden + torch.cat((xx, xx), dim=-1)
            else:
                xx = hidden
            
        return xx
                    
        
    def check(self):
        print("*" * 40)
        for tmp_i in range(1, len(self.sizes)):
            print(self.weights["w_{0}_{1}".format(tmp_i-1, tmp_i)].shape)
        print("*" * 40)
        
        print("*" * 40)
        for tmp_i in range(1, len(self.sizes)):
            print(self.bias["b_{0}_{1}".format(tmp_i-1, tmp_i)].shape)
        print("*" * 40)
        
        print("*" * 40)
        for tmp_i in range(1, len(self.sizes)):
            print(self.resnet_dts["r_{0}_{1}".format(tmp_i-1, tmp_i)].shape)
        print("*" * 40)