06.almtp, ... 28.almtp: are the standard MTP files with different number of paramters (levels 6 through 28)
The levels 2 and 4 are depreciated: their training may not be stable.
Some stuff should be manually edited, see the "<<<" comments

MTP
version = 1.1.0
potential_name = MTP1m
species_count = 3 <<< number of components
potential_tag = 
radial_basis_type = RBChebyshev
	min_dist = 2 <<< minimal reasonably possible distance, in angstrom (check the `mlp mindist` command)
	max_dist = 5 <<< maximal (cutoff) distance
	radial_basis_size = 8 <<< number of radial basis functions, typical values are (in order of frequency of use): 8,6,10,4,12
	radial_funcs_count = 2
alpha_moments_count = 18
alpha_index_basic_count = 11
alpha_index_basic = {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 2, 0, 0}, {0, 1, 1, 0}, {0, 1, 0, 1}, {0, 0, 2, 0}, {0, 0, 1, 1}, {0, 0, 0, 2}, {1, 0, 0, 0}}
alpha_index_times_count = 14
alpha_index_times = {{0, 0, 1, 11}, {1, 1, 1, 12}, {2, 2, 1, 12}, {3, 3, 1, 12}, {4, 4, 1, 13}, {5, 5, 2, 13}, {6, 6, 2, 13}, {7, 7, 1, 13}, {8, 8, 2, 13}, {9, 9, 1, 13}, {0, 10, 1, 14}, {0, 11, 1, 15}, {0, 12, 1, 16}, {0, 15, 1, 17}}
alpha_scalar_moments = 9
alpha_moment_mapping = {0, 10, 11, 12, 13, 14, 15, 16, 17}
