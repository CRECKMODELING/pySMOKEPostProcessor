// In the future think about removing this stuff from here
// include definition should be separated in each file

// Standard library
#include <sstream>
#include <vector>
#include <exception>
#include <stdexcept>
#include <string>
#include <iostream>
#include <numeric>

// OpenSMOKEpp library
#include "OpenSMOKE_Definitions.h"
#include "kernel/thermo/Thermodynamics_CHEMKIN"
#include "kernel/transport/Transport_CHEMKIN"
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"
#include "maps/FluxAnalysisMap.h"
#include "math/PhysicalConstants.h"
#include "kernel/kinetics/KineticsUtilityFunctions.h"

// Surface phase additions
#include "kernel/kinetics/ReactionPolicy_Surface_CHEMKIN.h"
#include "maps/ThermodynamicsMap_Surface_CHEMKIN.h"
#include "maps/KineticsMap_Surface_CHEMKIN.h"

// Boost library
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/program_options.hpp>

// Eigen
#include <Eigen/Dense>

// pyBIND11
#pragma once

#ifdef BUILD_PYTHON_BINDINGS
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
#endif

