#ifndef STRUCTURE_H
#define STRUCTURE_H


namespace matersdk {

template <typename CoordType>
class Structure {
    Structure(int num_atoms);
    Structure(int num_atoms, CoordType **basis_vectors, int *atomic_numbers, CoordType **coords);
    

private:
    CoordType **basis_vectors;
    int *atomic_numbers;
    CoordType **cart_coords;
}; // class: Structure

} // namespace: matersdk


#endif