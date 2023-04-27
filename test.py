import numpy as np
import multiprocessing as mp

# define the save function to be run in parallel
def save_array(filename):
    array = np.array([1, 2, 3])
    np.save(filename, array)


if __name__ == '__main__':
    # create a list of (array, filename) tuples
    filenames = [(f"array_{i}.npy", ) for i in range(5)]

    # create a Pool with 4 worker processes
    with mp.Pool(processes=4) as pool:
        # map the save function to the data and execute it in parallel
        pool.starmap(save_array, filenames)