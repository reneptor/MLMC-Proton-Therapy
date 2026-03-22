#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   // std::vector, std::string
#include <pybind11/numpy.h>
#include "simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(mlmcprotons, m) {

    py::class_<MLMCprotons<std::mt19937>>(m, "MLMCprotons")
        .def(py::init<unsigned int, const std::string&>(),
             py::arg("seed"), py::arg("tablePath"))

        .def("loadPhantom", &MLMCprotons<std::mt19937>::loadPhantom,
             py::arg("phantomPath"), py::arg("n_x"), py::arg("n_y"), py::arg("n_z"),
             "Load a phantom volume from file with given dimensions.")

        .def("addTreatmentPlan", &MLMCprotons<std::mt19937>::addTreatmentPlan,
             py::arg("level"),
             "Add a treatment plan for a given level.")

        .def("addPencilBeam", &MLMCprotons<std::mt19937>::addPencilBeam,
             py::arg("level"), py::arg("nPrimShare"), py::arg("initialStep"),
             py::arg("entranceDir") = dc::ENTRANCE_DIR, py::arg("beamWidth") = dc::BEAM_WIDTH,
             py::arg("dir_x") = dc::DIR_X, py::arg("dir_y") = dc::DIR_Y, py::arg("dir_z") = dc::DIR_Z,
             py::arg("x_0") = dc::ENTRY_X, py::arg("y_0") = dc::ENTRY_Y, py::arg("z_0") = dc::ENTRY_Z,
             py::arg("E") = dc::BEAM_ENERGY, py::arg("spread_E") = dc::BEAM_ENERGY_SPREAD,
             py::arg("alpha") = dc::ALPHA,
             "Add a pencil beam to a treatment plan at a given level.")

        .def("simulateTreatmentPlan", &MLMCprotons<std::mt19937>::simulateTreatmentPlan,
             py::arg("level"), py::arg("numPrimaries"), py::arg("numThreads") = MLMCprotons<std::mt19937>::maxThreads,
             "Simulate a treatment plan for a given level, optionally multithreaded.")

        .def("renderCombinedScoringGrid", &MLMCprotons<std::mt19937>::renderCombinedScoringGrid,
             "Combine all treatment plan grids into the combined scoring grid.")

        .def("yieldDoseAtLevel", &MLMCprotons<std::mt19937>::yieldDoseAtLevel,
             py::arg("level"),
             "Return a ScoringGrid (dose) at the specified level.")

        .def("yieldDoseCombined", &MLMCprotons<std::mt19937>::yieldDoseCombined,
             "Return the combined dose ScoringGrid.");

    py::class_<ScoringGrid>(m, "Dose")
        .def(py::init<unsigned int, size_t, size_t, size_t>(),
             py::arg("level"), py::arg("n_x"), py::arg("n_y"), py::arg("n_z"))

        .def_readwrite("level", &ScoringGrid::level)
        .def_readwrite("n_x", &ScoringGrid::n_x)
        .def_readwrite("n_y", &ScoringGrid::n_y)
        .def_readwrite("n_z", &ScoringGrid::n_z)

        .def_property_readonly("grid", [](ScoringGrid &self) {
            if (!self.grid)
                throw std::runtime_error("Grid already released or uninitialized");

            float* data = self.grid;
            self.grid = nullptr;  // avoid double-free

            ssize_t n_x = static_cast<ssize_t>(self.n_x);
            ssize_t n_y = static_cast<ssize_t>(self.n_y);
            ssize_t n_z = static_cast<ssize_t>(self.n_z);

            // Capsule to handle memory cleanup
            py::capsule free_when_done(data, [](void* ptr) {
                float* f = reinterpret_cast<float*>(ptr);
                delete[] f;
            });

            return py::array_t<float>(
                {n_x, n_y, n_z},
                {sizeof(float) * n_y * n_z, sizeof(float) * n_z, sizeof(float)},
                data,
                free_when_done
            );
        });
}
