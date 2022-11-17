from pymatgen.io.abinit.abiobjects import KSampling
from pymatgen.io.abinit.abiobjects import KSamplingModes



class KSampler(KSampling):
    '''
    self.__init__()

    Parameters
    ----------
        1. mode: `KSamplingModes.monkhorst`. 
            - Mode for generating k-points. Use one of the `KSamplingModes` enum types.
        2. num_kpts: `0`
            - Number of kpoints if mode is “automatic” Number of division for the sampling 
            of the smallest segment if mode is “path”. Not used for the other modes.
        3. kpts: ``
        4. kpt_shifts:
        5. kpts_weights:
        6. use_symmetries:
        7. use_time_reversal:
        8. chksymbreak:
        9. comment: 
    '''
    pass
