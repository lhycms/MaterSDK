import numpy as np


class TildeRPairNormalizer(object):
    def __init__(self, tildeRs_array:np.ndarray):
        # shape: (num_frames, num_centers, max_num_nbrs, 4)     e.g. (48, 26, 4)
        #   ->
        # shape: (num_frames * num_centers * max_num_nbrs, 4)   e.g. (1248, 4)
        self.tildeRs_array = tildeRs_array.reshape(-1, 4)
        # shape = (1, 4)
        self.davg, self.dstd = self.calc_stats()
    
    
    def calc_stats(self):
        '''
        Description
        -----------
            1. 计算 DeepPot-SE 中 TildeR 的平均值(`avg`)和方差(`std`)
        '''
        ### Step 1. 分别获取径向信息(`info_radius`)和角度信息(`info_angles`)
        # shape: (num_frames * num_centers * max_num_nbrs, 1)
        info_radius = self.tildeRs_array[:, 0].reshape(-1, 1)
        # shape: (num_frames * num_centers * max_num_nbrs, 3)
        info_angles = self.tildeRs_array[:, 1:].reshape(-1, 3)
        
        ### Step 2. 分别获取径向信息(`info_radius`)和角度信息(`info_angles`)的一些量:
        #     1) sum: 求和
        #     2) sum^2: 平方和 -- 先平方后求和
        #     3) total_num_pairs: num_centers * max_num_nbrs
        ### Step 2.1. sum: 求和
        sum_info_radius = np.sum(info_radius)
        sum_info_angles = np.sum(info_angles) / 3.0
        
        ### Step 2.2. sum^2: 平方和 -- 先平方后求和
        sum2_info_radius = np.sum(
                    np.multiply(info_radius, info_radius)
        )
        sum2_info_angles = np.sum(
                    np.multiply(info_angles, info_angles)
        ) / 3.0
        
        ### Step 2.3. total_num_pairs.shape: num_centers * max_num_nbrs
        total_num_pairs = info_radius.shape[0]
        
        
        ### Step 3. 计算平均值 -- davg_unit
        davg_unit = [sum_info_radius / (total_num_pairs + 1e-15), 0, 0, 0]
        # shape = (1, 4)
        davg_unit = np.array(davg_unit).reshape(1, 4)
        
        ### Step 4. 计算方差 -- dstd_unit
        dstd_unit = [
            self._calc_std(sum2_value=sum2_info_radius, sum_value=sum_info_radius, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs)
        ]
        # shape = (1, 4)
        dstd_unit = np.array(dstd_unit).reshape(1, 4)
        
        return davg_unit, dstd_unit
        
    
    def _calc_std(self, sum2_value:float, sum_value:float, N:int):
        '''
        Description
        -----------
            1. 计算标准差
        
        Parameters
        ----------
            1. sum2_value: float
                - sum2_value = \sum_i^N{x_i^2}，先平方后求和
            2. sum_value: float
                - sum_value  = \sum_i^N{x_i}
        '''
        if (N == 2):
            return 1e-2
        std = np.sqrt(
                sum2_value / N - np.multiply(sum_value/N, sum_value/N)
        )
        if np.abs(std) < 1e-2:
            std = 1e-2
        return std
        
    
    def normalize(self, tildeRs_array:np.ndarray):
        '''
        Parameters
        ----------
            1. tildeRs_array: np.ndarray
                - shape = (num_frames, num_centers, max_num_nbrs, 4)
        '''
        davg = self.davg.reshape(1, 1, 1, 4)
        dstd = self.dstd.reshape(1, 1, 1, 4)
        result = (tildeRs_array - davg) / dstd
        
        return result
        