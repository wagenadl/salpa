// salpapy.cpp

#include "LocalFit.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(salpa_cppcore, m) {
  m.doc() = "SALPA plugin";
  py::class_<LocalFit>(m, "Salpa")
    //    .def(py::init<raw_t const *, raw_t *, timeref_t, raw_t, int>())
    .def("reset", &LocalFit::reset)
    .def("set_t_blankdepeg", &LocalFit::set_t_blankdepeg)
    .def("set_t_ahead", &LocalFit::set_t_ahead)
    .def("set_t_chi2", &LocalFit::set_t_chi2)
    .def("setrail", &LocalFit::setrail)
    .def("process", &LocalFit::process)
    .def("forcepeg", &LocalFit::forcepeg);
  m.def("isvalid", []() { return LocalFit::isValid(); });
  m.def("csalpa", [](py::array_t<float_t, py::array::c_style> in,
                     py::array_t<float, py::array::c_style> out,
                     float thresh, int tau) {
    auto in_ = in.request();    
    auto out_ = out.request();    
    return new LocalFit(static_cast<float *>(in_.ptr),
                        static_cast<float *>(out_.ptr),
                        in_.shape[0],
                        thresh,
                        tau); });
}


