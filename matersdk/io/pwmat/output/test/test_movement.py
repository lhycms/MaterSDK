import unittest

# python3 -m matersdk.io.pwmat.output.test.test_movement
from ..movement import Movement



class MovementTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        movement = Movement(movement_path=movement_path)
        idx_frame = 300 # 帧数从 0 开始计数
        
        # 1. get_chunksize()
        print()
        print("Step 1. The chunksize of each frame:", end="\t")
        print(movement.get_chunksize())
        #print(movement.get_chunkslice())
        
        # 2. Get chunk of specific frame structure
        print()
        print("Step 2. String for structure:")
        atom_config_string = movement._get_frame_str(idx_frame=idx_frame)
        print(atom_config_string)
        
        # 3. 
        #print()
        #print("Step 3. Structure from {0}th frame:".format(idx_frame))
        #print(movement.get_frame_structure(idx_frame=idx_frame))
        
        '''
        # 4. 
        print()
        print("Step 4. Virial tensor of {0}th frame:".format(idx_frame)) 
        frame_structure = movement.get_frame_structure(idx_frame=idx_frame)
        print(frame_structure.virial_tensor)   
        
        # 5.
        print()
        print("Step 5. Energy of {0}th frame:".format(idx_frame))
        energy_tot, energy_p, energy_k = movement.get_frame_energy(idx_frame=idx_frame)
        print( "\t1. Total energy: {0} eV".format(energy_tot) )
        print( "\t2. Potential energy: {0} eV".format(energy_p) )
        print( "\t3. Kenitic energy: {0} eV".format(energy_k) )
        
        
        # 6. 
        print()
        print("Step 6. Virial tensor of {0}th frame:".format(idx_frame))
        virial_tensor = movement.get_frame_virial(idx_frame=idx_frame)
        print(virial_tensor)
        

        # 7. 
        print()
        print("Step 7. Volume of {0}th frame:`".format(idx_frame), end="\t")
        volume = movement.get_frame_volume(idx_frame=idx_frame)
        print(volume)
        
        
        # 8.
        print()
        print("Step 8. Force of {0}th frame:`".format(idx_frame), end="\t")
        forces_array = movement.get_frame_structure(idx_frame=idx_frame).get_atomic_force()
        print(forces_array)
        
        # 9.
        print()
        print("Step 9. Atomic energy of {0}th frame:`".format(idx_frame), end="\t")
        forces_array = movement.get_frame_structure(idx_frame=idx_frame).get_atomic_energy()
        print(forces_array)
        '''
    
        
        
if __name__ == "__main__":
    unittest.main()