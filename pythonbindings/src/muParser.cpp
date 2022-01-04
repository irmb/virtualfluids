#include <pybind11/pybind11.h>
#include <muParser.h>

namespace py = pybind11;

PYBIND11_MODULE(pymuparser, m) {
    py::class_<mu::ParserBase>(m, "_ParserBase");

    py::class_<mu::Parser, mu::ParserBase>(m, "Parser")
            .def(py::init())
            .def_property("expression", &mu::Parser::GetExpr, &mu::Parser::SetExpr)
            .def("set_decimal_separator", &mu::Parser::SetDecSep)
            .def("define_constant", &mu::Parser::DefineConst)
            .def("define_variable", &mu::Parser::DefineVar);

}