import h5py
import torch
import numpy as np


class FFExtractor(object):
    def __init__(self, 
                pt_path:str,
                hdf5_path:str,
                num_embedding_nets:int,
                num_embedding_layers:int,
                num_fitting_nets:int,
                num_fitting_layers:int):
        '''
        Parameters
        ----------
            1. pt_path: str
                - The path of checkpoint file.
            2. hdf5_path: str
                - The path of checkpoint file.
            3. num_embedding_nets: int
                - 没有AtomEmbedding的情况下, `num_embedding_nets` = 元素种类^2
            4. num_embedding_layers: int
                - embedding 的层数
            5. num_fitting_nets: int
                - 没有AtomEmbedding的情况下, `num_fitting_nets` = 元素种类
            6. num_fitting_layers: int
                - fitting layer 的层数
        '''
        self.model = torch.load(f=pt_path, map_location=torch.device("cpu"))
        self.model_state_dict = self.model["state_dict"]
        self.hdf5_path = hdf5_path
        
        self.num_embedding_nets = num_embedding_nets
        self.num_embedding_layers = num_embedding_layers
        self.num_fitting_nets = num_fitting_nets
        self.num_fitting_layers = num_fitting_layers


    def save_hdf5_file(self):
        '''
        Description
        -----------
            1. 获取 force field 的参数：
                1. weights
                2. bias
                3. resnet block
        '''
        ### Step 1. Create a HDF5 file
        hdf5_file = h5py.File(self.hdf5_path, 'w')
        
        ### Step 2. Extract force field parameters
        ### Step 2.1. Embedding net (Note: Resnet block)
        for tmp_net_idx in range(self.num_embedding_nets):
            for tmp_layer_idx in range(self.num_fitting_layers):
                hdf5_file.create_dataset(
                    "embedding_net.{0}.weights.weight{1}".format(tmp_net_idx, tmp_layer_idx),
                    data=self.model_state_dict["embedding_net.{0}.weights.weight{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                )
                hdf5_file.create_dataset(
                    "embedding_net.{0}.bias.bias{1}".format(tmp_net_idx, tmp_layer_idx),
                    data=self.model_state_dict["embedding_net.{0}.bias.bias{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                )
                ### Note: Resnet block of Embedding Net
                try:
                    hdf5_file.create_dataset(
                        "embedding_net.{0}.resnet_dt.resnet_dt{1}".format(tmp_net_idx, tmp_layer_idx),
                        data=self.model_state_dict["embedding_net.{0}.resnet_dt.resnet_dt{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                )
                except KeyError:
                    pass
        
        ### Step 2.2. Fitting net (Note: Resnet block)
        for tmp_net_idx in range(self.num_fitting_nets):
            for tmp_layer_idx in range(self.num_fitting_layers):
                hdf5_file.create_dataset(
                    "fitting_net.{0}.weights.weight{1}".format(tmp_net_idx, tmp_layer_idx),
                    data=self.model_state_dict["embedding_net.{0}.weights.weight{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                )
                hdf5_file.create_dataset(
                    "fitting_net.{0}.bias.bais{1}".format(tmp_net_idx, tmp_layer_idx),
                    data=self.model_state_dict["fitting_net.{0}.bias.bias{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                )
                ### Note: Resnet block of Fitting Net
                try:
                    hdf5_file.create_dataset(
                        "fitting_net.{0}.resnet_dt.resnet_dt{1}".format(tmp_net_idx, tmp_layer_idx),
                        data=self.model_state_dict["fitting_net.{0}.resnet_dt.resnet_dt{1}".format(tmp_net_idx, tmp_layer_idx)].numpy()
                    )
                except KeyError:
                    pass
        
        
        ### Step 3. Close a HDF5 file
        hdf5_file.close()
        
        
    @staticmethod
    def read_hdf5_file(hdf5_path:str):
        pass
    
    
    @staticmethod
    def get_params():
        pass