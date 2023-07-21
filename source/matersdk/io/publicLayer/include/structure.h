#ifndef STRUCTURE_H
#define STRUCTURE_H


namespace matersdk {

template <typename CoordType>
class Structure {
    Structure(int num_atoms);
    Structure(int num_atoms, CoordType **basis_vectors, int *atomic_numbers, CoordType **coords);
    Structure(const Structure &rhs);
    ~Structure();
    //void make_supercell(int *scaling_matix);


private:
    int num_atoms;
    CoordType **basis_vectors;
    int *atomic_numbers;
    CoordType **cart_coords;
}; // class: Structure

} // namespace: matersdk


#endif