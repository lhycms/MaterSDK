#include <pybind11/pybind11.h>
#include "../include/se.h"


PYBIND11_MODULE(se, m) {
    py::class_<matersdk::deepPotSE::TildeR>()
        .def(py::init<int>())
}
