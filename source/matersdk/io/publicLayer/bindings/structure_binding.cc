#include <pybind11/pybind11.h>
#include "../include/structure.h"

namespace py = pybind11;


PYBIND11_MODULE(Structure, m) {
    py::class_<matersdk::Structure<double>>(m, "Structure")
        .def(py::init<>)
}