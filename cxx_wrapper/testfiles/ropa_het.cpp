#include "../PostProcessorWrapper_py.h"
#include "core/PostProcessorCore.h"
#include "maps/Maps_CHEMKIN"
#include "ropa/ROPA.h"

int main(int argc, char** argv) {
  const std::string mechanism_folder =
      std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Surface_Data/kinetics";
  const std::string output_folder =
      std::string(PYSMOKE_EXAMPLES_DATA_DIR) + "/Surface_Data/Output";

  PostProcessorCore profiles_db;
  bool dummy = profiles_db.ReadOutput(output_folder, true);
  dummy = profiles_db.LoadKinetics(mechanism_folder);

  const unsigned int N_rxns = 5;
  std::vector<unsigned int> rxns_indices(N_rxns);
  std::vector<double> rxns_coefficients(N_rxns);

  bool heterogeneous_reactions = false;

  ROPA ropa;
  ropa.SetResults(&profiles_db);
  ropa.SetROPAType("global");
  // ropa.SetSpecies("CH4");
  // ropa.RateOfProductionAnalysis(N_rxns,heterogeneous_reactions);

  // rxns_indices = ropa.reactions();
  // rxns_coefficients = ropa.coefficients();

  // std::cout << "Species CH4" << std::endl;
  // std::cout << "Homogeneous reactions " << std::endl;
  // for (unsigned int i = 0; i<N_rxns; i++){
  //     std::cout << "Reaction index: " << rxns_indices[i] << std::endl;
  //     std::cout << "Reaction coeff: " << rxns_coefficients[i] << std::endl;
  // }

  // Setting heterogeneous reactions
  heterogeneous_reactions = true;
  // ropa.RateOfProductionAnalysis(N_rxns,heterogeneous_reactions);

  // rxns_indices = ropa.reactions();
  // rxns_coefficients = ropa.coefficients();

  // std::cout << "Heterogeneous reactions " << std::endl;
  // for (unsigned int i = 0; i<N_rxns; i++){
  //     std::cout << "Reaction index: " << rxns_indices[i] << std::endl;
  //     std::cout << "Reaction coeff: " << rxns_coefficients[i] << std::endl;
  // }

  ////////////////////////////// CH3(S)
  // std::cout << std::endl;
  // ropa.SetSpecies("CH3(S)");

  // std::cout << "Species CH3(S)" << std::endl;
  // ropa.RateOfProductionAnalysis(N_rxns,heterogeneous_reactions);

  // rxns_indices = ropa.reactions();
  // rxns_coefficients = ropa.coefficients();

  // std::cout << "Heterogeneous reactions " << std::endl;
  // for (unsigned int i = 0; i<N_rxns; i++){
  //     std::cout << "Reaction index: " << rxns_indices[i] << std::endl;
  //     std::cout << "Reaction coeff: " << rxns_coefficients[i] << std::endl;
  // }

  //////////////////////// Species C(B)
  std::string sp_name = "CH3(S)";

  std::cout << std::endl;
  ropa.SetSpecies(sp_name);
  std::cout << "Species " << sp_name << std::endl;
  ropa.RateOfProductionAnalysis(N_rxns, heterogeneous_reactions);

  rxns_indices = ropa.reactions();
  rxns_coefficients = ropa.coefficients();

  std::cout << "Heterogeneous reactions " << std::endl;
  for (unsigned int i = 0; i < N_rxns; i++) {
    std::cout << "Reaction index: " << rxns_indices[i] << std::endl;
    std::cout << "Reaction coeff: " << rxns_coefficients[i] << std::endl;
  }

  return 0;
}