import torch
import torch.nn as nn

from typing import List


class DpEmbeddingNet(nn.Module):
    def __init__(
                self,
                embedding_sizes:List[int],
                activation_fn=nn.Tanh,
                ):
        super(DpEmbeddingNet, self).__init__()
        self.sizes = [1] + embedding_sizes
        self.activation_fn = activation_fn()
        self.weights_dict = nn.ParameterDict()
        self.bias_dict = nn.ParameterDict()
        
        ### 初始化 weights 和 bias
        for tmp_i in range(1, len(self.sizes)):
            ### weights
            tmp_weights = torch.Tensor(self.sizes[tmp_i-1], self.sizes[tmp_i])
            nn.init.normal_(tmp_weights, mean=0, std=0.01)
            self.weights_dict["w_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_weights)
            
            ### bias
            tmp_bias = torch.Tensor(self.sizes[tmp_i])
            nn.init.normal_(tmp_bias, mean=0, std=0.01)
            self.bias_dict["b_{0}_{1}".format(tmp_i-1, tmp_i)] = nn.Parameter(data=tmp_bias)
        
    
    def forward(self, tilde_r:torch.Tensor):
        x = tilde_r[:, :, :, 0].unsqueeze(dim=-1)
        for tmp_i in range(1, len(self.sizes)):
            x = torch.matmul(x, self.weights_dict["w_{0}_{1}".format(tmp_i-1, tmp_i)])
            x = x + self.bias_dict["b_{0}_{1}".format(tmp_i-1, tmp_i)].view(1, -1)
            x = self.activation_fn(x)
        g = x
        
        D = torch.matmul(g.transpose(-2, -1), tilde_r)
        D = torch.matmul(D, tilde_r.transpose(-2, -1))
        D = torch.matmul(D, g[:, :, :, :4])
        
        return D
    
    
    def check(self):
        for tmp_i in range(1, len(self.sizes)):
            print('*' * 40)
            print("w_{0}_{1} = ".format(tmp_i-1, tmp_i), self.weights_dict["w_{0}_{1}".format(tmp_i-1, tmp_i)].size())
            print("b_{0}_{1} = ".format(tmp_i-1, tmp_i), self.bias_dict["b_{0}_{1}".format(tmp_i-1, tmp_i)].size())
            print('*' * 40)
