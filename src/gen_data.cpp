/*-------------------------------------------------------------------
 * Brief description of this file: 
 * testing the model output of the reduced 5 rxn mechanism 
 *-----------------------------------------------------------------*/

#include <random>
#include "model.h"
#include "reaction_info.h"
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/cea_mixture.h>
#include <antioch/cea_evaluator.h>
#include <antioch/reaction_set.h>
#include <antioch/xml_parser.h>
#include <antioch/ascii_parser.h>
#include <antioch/read_reaction_set_data.h>
#include <antioch/cea_mixture_ascii_parsing.h>

int main()
{

  unsigned int n_species = 9;
  const int dim = n_species + 1;
  std::vector<double> initialValues(dim,0.); 
  std::vector<double> phiPoints(1,0.); 
  phiPoints[0] = 1.0; //phiPoints[1] = 1.0; phiPoints[2] = 1.1;
  std::vector<double> tempPoints(1,0.); 
  //double temp1 = 1300.;
  tempPoints[0] = 1200.;
  /* tempPoints[1] = 1300.; */
  /* tempPoints[2] = 1400.; */
  //double timePoint = 2.0e-5;

  std::vector<double> timePoints(10,0.);
  for (unsigned int i = 0; i < timePoints.size(); i++){
    timePoints[i] = (i+1) * 2.e-5/2;
    /* timePoints[i] = (i) * 2.5e-5 + 3.0e-5; */
    /* timePoints[i] = (i+1) * 5.0e-6; */
  }
  
  //open file now to write data
  std::ofstream datafile;
  datafile.open ("/home/rebecca/repos/combustion/time-series-data/datafile.txt");
  
  for (unsigned int l = 0; l < phiPoints.size(); l++){
   for (unsigned int ll = 0; ll < tempPoints.size(); ll++){
    initialValues[0] = 2.0*phiPoints[l]; //stoichiometric
    initialValues[1] = 1.0; //O2
    initialValues[2] = 0.0; //H
    initialValues[3] = 0.0; //O
    initialValues[4] = 0.0; //OH
    initialValues[5] = 0.0; //HO2
    initialValues[6] = 0.0; //H2O
    initialValues[7] = 0.0; //H2O2
    initialValues[8] = 3.78; //N2
    initialValues[9] = tempPoints[ll]; //temp
    
    //set up antioch
    std::string input_name = "/home/rebecca/repos/combustion/xml-inputs/h2o2-21.xml";
  
    std::vector<std::string> species_str_list;
    species_str_list.reserve(n_species);
    species_str_list.push_back("H2");
    species_str_list.push_back("O2");
    species_str_list.push_back("H");
    species_str_list.push_back("O");
    species_str_list.push_back("OH");
    species_str_list.push_back("HO2");
    species_str_list.push_back("H2O");
    species_str_list.push_back("H2O2");
    species_str_list.push_back("N2");
  
    //return 10 time points of all species and temp
    std::vector<double> returnValues(dim*timePoints.size(),0.);
  
    //const double Tref = 1.0;//
    Antioch::ChemicalMixture<double> chem_mixture( species_str_list );
    Antioch::CEAThermoMixture<double> cea_mixture( chem_mixture );
    Antioch::ReactionSet    <double> reaction_set( chem_mixture );
  
    Antioch::read_reaction_set_data_xml<double>( input_name, true, reaction_set );
    Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
  
    //double R = Antioch::Constants::R_universal<double>();
  
    Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0);
    Antioch::CEAEvaluator<double> thermo(cea_mixture);
  
    std::vector<double> scales(1, 0.); //just need a placeholder here

    reaction_info rxnMain(&n_species, &chem_mixture, &reaction_set, &thermo, &kinetics, scales);
  
    hydrogenComputeModel(initialValues,timePoints,&rxnMain,returnValues);
   
    //create measurement error
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,std::sqrt(.001));
    //TODO increase this appropriately for temperature measurements

    //right now not adding any error
    for (unsigned int i = 0; i < timePoints.size(); i++){
      datafile << phiPoints[l] << " " << tempPoints[ll] << " " << timePoints[i] << " ";
      for (unsigned int j = 0; j < dim; j++){
          double error = distribution(generator);
          double measurement = returnValues[dim*i+j];//add error:// + error;
          if (measurement>0){datafile << measurement << " ";}
          else {datafile << 0.0 << " ";}
      } datafile << "\n";
    }
 
   }
  }
  return 0;
}
