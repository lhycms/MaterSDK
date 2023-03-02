from typing import List


class DosInput(object):
    '''
    Description
    -----------
        1. 使用 `plot_DOS.interp.x` 的时候，需要有 `DOS.input` 作为输入
        2. `DOS.input` 的内容:
            ```
            0   # setting 0 -> for all atoms; setting 1 -> for partial atoms. 
            0   # setting 1 -> using interpolation; setting 0 -> using Gaussian broadening
            0.05    4000    #  energy smearing, in eV;  number of energy grid points, default is 4000
            8   8   8   # interpolation grid
            ```
    '''
    def __init__(
            self,
            mark_atoms:int,
            mark_method:int,
            ismear:float, num_energies:int,
            grid_interp:List[float],
            ):
        self.mark_atoms = mark_atoms
        self.mark_method = mark_method
        self.ismear = ismear
        self.num_energies = num_energies
        self.grid_interp = grid_interp
        
    
    def to(self, output_path:str):
        '''
        Desscription
        ------------
            1. 生成 `DOS.input`
        '''
        with open(output_path, "w") as f:
            f.write("{0}\n".format(self.mark_atoms))
            f.write("{0}\n".format(self.mark_method))
            f.write("{0:<.5f}  {1:<10d}\n".format(self.ismear, self.num_energies))
            f.write("  ".join([str(value) for value in self.grid_interp]))
            