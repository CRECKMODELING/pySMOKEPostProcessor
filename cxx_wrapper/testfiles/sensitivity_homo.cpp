#include "../PostProcessorWrapper_py.h"
#include "core/PostProcessorCore.h"
#include "sensitivity/SensitivityReader.h"
//#include "ROPA.h"
#include "maps/Maps_CHEMKIN"

int main(int argc, char** argv)
{
    const std::string mechanism_folder = std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Sensitivity/kinetics";
    const std::string output_folder = std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Sensitivity/Output";

    PostProcessorCore profiles_db;
    bool dummy = profiles_db.ReadOutput(output_folder, false);
    dummy = profiles_db.LoadKinetics(mechanism_folder);

    SensitivityReader sensi;
    sensi.SetResults(&profiles_db);
    sensi.SetSensitivityType("global");
    sensi.SetOrderingType("peak-values");
    sensi.SetNormalizationType("max-value");
    sensi.SetTarget("H2");
    sensi.Prepare();
    sensi.ReadSensitivityCoefficients();
    sensi.Sensitivity_Analysis(5);
    std::vector<unsigned int> reactions = sensi.reactions();
    std::vector<double> sensicoeffs = sensi.sensitivityCoefficients();
    for (unsigned int i = 0; i<sensicoeffs.size(); i++)
        std::cout << sensicoeffs[i] << std::endl;

    return 0;
}