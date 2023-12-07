#ifndef PYMATGEN_UTILS_LOAD_FILE_H

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>

namespace PymatgenUtils {
class Structure {
public:
    Structure();

    int from_file(const char* file_path);

    //Structure(const Structure& rhs);

    //Structure& operator=(const Structure& rhs);

    //~Structure();

    //void find_lattice(double** lattice);

    //void find_atomic_numbers(int* atomic_numbers);

    //void find_cart_coords(double** cart_coords);

private:
    //PyObject* structure_py = nullptr;
};  // class : StructureLoader


}   // namespace : PymatgenUtils

#endif