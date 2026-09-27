#include "../PostProcessorWrapper_py.h"
#include "../source/core/PostProcessorCore.h"
#include "../source/sensitivity/SensitivityCalculator.h"
// #include "ROPA.h"
#include "maps/Maps_CHEMKIN"

int main(int argc, char** argv) {
  const std::string mechanism_folder =
      std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Surface_Data/kinetics";
  const std::string output_folder =
      std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Surface_Data/Output";

  PostProcessorCore profiles_db;
  bool dummy = profiles_db.ReadOutput(output_folder, true);
  dummy = profiles_db.LoadKinetics(mechanism_folder);
  SensitivityCalculator sensi;
  const unsigned int NR = 5;

  sensi.SetResults(&profiles_db);
  sensi.SetSensitivityType("global");
  sensi.SetOrderingType("peak-values");
  sensi.SetNormalizationType("max-value");
  sensi.SetTarget("CH4");
  // Until here, everything is the same for sensitivity homogeneous and heterogeneous

  // Homogeneous sensitivity analysis
  sensi.Prepare(Phase::Gas);
  sensi.ReadSensitivityCoefficients();
  sensi.Sensitivity_Analysis(NR);
  std::vector<unsigned int> reactions = sensi.reactions();
  std::vector<double> sensicoeffs = sensi.sensitivityCoefficients();

  std::cout << "Homogeneous Sensitivity coefficients" << std::endl;
  for (unsigned int i = 0; i < sensicoeffs.size(); i++)
    std::cout << sensicoeffs[i] << std::endl;

  std::cout << "Indices" << std::endl;
  for (unsigned int i = 0; i < reactions.size(); i++)
    std::cout << reactions[i] << std::endl;

  // Heterogeneous sensitivity analysis
  sensi.Prepare(Phase::Surface);
  sensi.ReadSensitivityCoefficients();
  sensi.Sensitivity_Analysis(NR);
  std::vector<unsigned int> reactions_het = sensi.reactions();
  std::vector<double> sensicoeffs_het = sensi.sensitivityCoefficients();

  std::cout << "\nHeterogeneous Sensitivity coefficients" << std::endl;
  for (unsigned int i = 0; i < sensicoeffs_het.size(); i++)
    std::cout << sensicoeffs_het[i] << std::endl;

  std::cout << "Indices" << std::endl;
  for (unsigned int i = 0; i < reactions.size(); i++)
    std::cout << reactions_het[i] << std::endl;

  return 0;
}