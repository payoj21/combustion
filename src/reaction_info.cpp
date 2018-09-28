#include "reaction_info.h"

//Constructor
reaction_info::reaction_info(
    const unsigned int * n_species,
    Antioch::ChemicalMixture<double> * chem_mixture,
    Antioch::ReactionSet<double> * reaction_set, 
    Antioch::CEAEvaluator<double> * thermo,
    Antioch::KineticsEvaluator<double> * kinetics,
    std::vector<double> & scales)
:
  N_species(n_species),
  Chem_mixture(chem_mixture),
  Reaction_set(reaction_set),
  Thermo(thermo),
  Kinetics(kinetics),
  Scales(scales)
{
}

//Destructor
reaction_info::~reaction_info()
{
}
