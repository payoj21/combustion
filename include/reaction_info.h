#ifndef REACTION_INFO_H
#define REACTION_INFO_H

//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/reaction_set.h>
#include <antioch/cea_evaluator.h>
#include <antioch/kinetics_evaluator.h>
//c++
#include <vector>

// define struct that holds all reaction info, except params
struct reaction_info { reaction_info(
  const unsigned int * n_species,
  Antioch::ChemicalMixture<double> * chem_mixture,
  Antioch::ReactionSet<double> * reaction_set, 
  Antioch::CEAEvaluator<double> * thermo,
  Antioch::KineticsEvaluator<double> * kinetics,
  std::vector<double> & scales);
 ~reaction_info();

  const unsigned int * N_species;
  Antioch::ChemicalMixture<double> * Chem_mixture;
  Antioch::ReactionSet<double> * Reaction_set;
  Antioch::CEAEvaluator<double> * Thermo;
  Antioch::KineticsEvaluator<double> * Kinetics;
  std::vector<double> & Scales;
};
#endif
