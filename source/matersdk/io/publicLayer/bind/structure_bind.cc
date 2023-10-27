#include <pybind11/pybind11.h>
#include "../include/structure.h"



PYBIND11_MODULE(structure, m) {
    pybind11::class_<matersdk::Structure<double>>(m, "CStructure")
        //.def(pybind11::init<>());
        .def(
            pybind11::init<int, double**, int*, double**, bool>(), 
            pybind11::arg("num_atoms"), 
            pybind11::arg("basis_vectors"), 
            pybind11::arg("atomic_numbers"), 
            pybind11::arg("coords"), 
            pybind11::arg("is_cart_coords")=true
        );
        
    
}