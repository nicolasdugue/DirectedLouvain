/* Based on the following sources:
 * - https://ekbanaml.github.io/old_site_data/cpp-references/pybind-reference/
 * - https://pybind11.readthedocs.io/en/latest/basics.html#default-args
 */
#include "../include/community.hpp"
#include "../extern/pybind11/include/pybind11/pybind11.h"
#include "../extern/pybind11/include/pybind11/operators.h"
#include "../extern/pybind11/include/pybind11/iostream.h"
#include "../extern/pybind11/include/pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(directedlouvain, dl) {
    py::class_<Community>(dl, "Community")
            .def(py::init<const string &, bool, const double, const double, bool, bool>(),
                 py::arg("filename"), py::arg("weighted") = false, py::arg("precision") = 0.0001, py::arg("gamma") = 1,
                 py::arg("reproducibility") = false, py::arg("renumbering") = true)
            .def("run", &Community::run,
                 py::arg("verbose") = false, py::arg("display_level") = -1, py::arg("filename_part") = "",
                 py::call_guard<py::scoped_ostream_redirect,
                 py::scoped_estream_redirect>())
            .def("modularity", &Community::modularity)
            .def("last_level", &Community::get_last_level)
            .def("get_level", &Community::get_level,
                 py::arg("level") = 0)
            .def("print_level", &Community::print_level,
                 py::arg("level") = 0,
                 py::call_guard<py::scoped_ostream_redirect,
                 py::scoped_estream_redirect>());
}
