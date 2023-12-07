#include <iostream>

#include "../include/load_file.h"


namespace PymatgenUtils {
/**
 * @brief Default Construct : do nothing
 * 
 */
Structure::Structure() {
}


int Structure::from_file(const char* file_path) 
{   
    Py_Initialize();
    import_array();
    // Step 1. Import `pymatgen.core` module, get `Structure.from_file()` method
    PyObject* pymatgen_core_str_py = PyUnicode_FromString("pymatgen.core");
    PyObject* pymatgen_core_module_py = PyImport_Import(pymatgen_core_str_py);
    Py_DECREF(pymatgen_core_str_py);
    PyObject* structure_class_py = PyObject_GetAttrString(pymatgen_core_module_py, "Structure");
    Py_DECREF(pymatgen_core_module_py);
    PyObject* from_file_method_py = PyObject_GetAttrString(structure_class_py, "from_file");
    Py_DECREF(structure_class_py);

    // Step 2. run `Structure.from_file()`
    PyObject* args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(file_path));
    PyObject* structure_py = PyObject_CallObject(from_file_method_py, args);
    Py_DECREF(from_file_method_py);
    Py_DECREF(args);
    if (!structure_py) {
        std::cerr << "Can't read structure file!!!";
    }

    // Step 3. 
    //this->structure_py = structure_py;
    Py_DECREF(structure_py);

    Py_Finalize();
    
    return 0;
}


}   // namespace : PymatgenUtils