 /*------------------------------------------------------------------
 * This file is divided in two parts:
 * - the first one handles the statistical inverse problem (SIP) for estimating
 *   the inadequacy parameters
 * - the second part handles the statistical forward problem (SFP) for
 *   predicting the QoI
 *
 * The SIP definition requires a user defined likelihood function; refer to files 
 * 'likelihood.h' and 'likelihood.cpp'. The SFP definition requires a user 
 * defined qoi function; refer to files 'qoi.h' and 'qoi.cpp'.
 *-----------------------------------------------------------------*/

#include "compute.h"
#include "likelihood.h"
#include "qoi.h"
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
//queso
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/UniformVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <sys/time.h>
#include <cmath>

void computeReactionRates(const QUESO::FullEnvironment& env) {

  struct timeval timevalNow;
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\nBeginning run of 'H2O2 inverse problem' at "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  //================================================================
  // Statistical inverse problem (SIP)
  //================================================================
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> Inadequacy parameter estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  const unsigned int n_species = 9;                 //the number of species included in the model
  const unsigned int dim = n_species + 1;           //the total dimension for the ODE solve, species + temperature
  const unsigned int n_reactions = 5;               //the number of reactions in reduced model 
  const unsigned int n_params = 3 * n_reactions;    //the number of parameters in the reduced model (A, n, E), k = A T^n exp(-E/RT)

  //------------------------------------------------------
  // SIP Step 0 of 6: Read in the data
  //------------------------------------------------------
  double tmpPhis;
  double tmpTICs;
  double tmpTimes;
  double tmpH2;
  double tmpO2;
  double tmpH;
  double tmpO;
  double tmpOH;
  double tmpHO2;
  double tmpH2O;
  double tmpH2O2;
  double tmpN2;
  double tmpTemps;
  int numLines = 0;
  int numData = 0;

  FILE *dataFile;
  dataFile = fopen("/home/rebecca/repos/combustion/time-series-data/datafile.txt","r");

  std::vector<double> phis;//n_phis * n_tICs * n_times, 0.);
  std::vector<double> tICs;//n_phis * n_tICs * n_times, 0.);
  std::vector<double> times;//n_phis * n_tICs * n_times, 0.);
  std::vector<double> concs;//n_phis * n_tICs * n_times * n_species, 0.);
  std::vector<double> temps;//n_phis * n_tICs * n_times, 0.);
  /* std::cout << "HELLO? " << std::endl; */
  while (fscanf(dataFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tmpPhis, &tmpTICs,
        &tmpTimes, &tmpH2, &tmpO2, &tmpH, &tmpO, &tmpOH, &tmpHO2, &tmpH2O, &tmpH2O2, &tmpN2, &tmpTemps) != EOF) {
    phis.push_back(tmpPhis);
    tICs.push_back(tmpTICs);
    times.push_back(tmpTimes);
    concs.push_back(tmpH2);
    concs.push_back(tmpO2);
    concs.push_back(tmpH);
    concs.push_back(tmpO);
    concs.push_back(tmpOH);
    concs.push_back(tmpHO2);
    concs.push_back(tmpH2O);
    concs.push_back(tmpH2O2);
    concs.push_back(tmpN2);
    temps.push_back(tmpTemps);
    numLines++;
    numData+=n_species;
  }

  fclose(dataFile);
  std::cout << "The number of data points is " << numData << std::endl;
  int i=0;
  while (times[i+1] > times[i]){
      i++;
  }
  const unsigned int n_times = i;
  std::cout << "The number of initial conditions is " << numLines/n_times << std::endl;
  std::cout << "The number of time points is " << n_times << std::endl;

  const unsigned int n_ICs = numLines/n_times;                     //the number of initial conditions (equivalence ratio)
  //set up antioch
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

  //const double Tref = 1.0;//
  std::cout<<species_str_list<<std::endl;
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );
  std::cout << " HI  \n" << std::endl;
  Antioch::CEAThermoMixture<double> cea_mixture( chem_mixture );
  Antioch::ReactionSet    <double> reaction_set( chem_mixture );

  std::string kinetics_input_name = "xml-inputs/h2o2-5.xml";
  Antioch::read_reaction_set_data_xml<double>( kinetics_input_name, true, reaction_set );
  Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );

  //double R = Antioch::Constants::R_universal<double>();

  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0);
  Antioch::CEAEvaluator<double> thermo(cea_mixture);

  //std::cout << "rxnMain.Reaction_set.reaction(1) is " << rxnMain.Reaction_set->reaction(1) << std::endl;

  //------------------------------------------------------
  // SIP Step 1 of 6: Instantiate the parameter space
  //------------------------------------------------------
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", n_params, NULL);
  
  //------------------------------------------------------
  // SIP Step 2 of 6: Instantiate the parameter domain
  //------------------------------------------------------
  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());
  
  //paramMinValues.cwSet(0.);
  //paramMaxValues.cwSet(10.);
  // For uniform priors:
  /*double min_temp[n_params] = {2.5, 3.5, 0.9, 5., 1., -1., 2., 1., -2., -1., 60., 20., 10., 0., 240.};
  double max_temp[n_params] = {4.5, 6.5, 1.3, 6., 2., 0.1, 3., 2.5, -1., 1.1, 80., 30., 20., 1., 270.};
  for (unsigned int i=0; i<n_params; i++){
    paramMinValues[i] = min_temp[i];
    paramMaxValues[i] = max_temp[i];
  }*/
 //For gaussian priors:
  // A (prefactor) have lognormal prior, n and E have gaussian priors
  // must take exponential of A params later
  for (unsigned int i=0; i<n_params; i++){
    paramMinValues[i] = -INFINITY;
    paramMaxValues[i] =  INFINITY;
  }
  //End Gaussian priors

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);
  
  //scaling factors used for parameters
  double scaling[n_params] = {
    1.0e10, 1.0e-2, 1.0e3, 1.0e7, 1.0e8,
    1.0e0,  1.0e0,  1.0e0, 1.0e0, 1.0e0, 
    1.0e3,  1.0e3,  1.0e3, 1.0e3, 1.0e3};
  std::vector<double> scales(n_params,0.0);
  for (unsigned int i = 0; i < n_params; i++) {scales[i] = scaling[i];}
  
  reaction_info rxnMain(&n_species, &chem_mixture, &reaction_set, &thermo, &kinetics, scales);

  //------------------------------------------------------
  // SIP Step 3 of 6: Instantiate the likelihood function 
  // object to be used by QUESO.
  //------------------------------------------------------

  //right now temperature is not being passed in, but it could be
  likelihoodRoutine_Data likelihoodRoutine_Data1(env, phis, tICs, times, concs, temps, &rxnMain);

  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj(
        "like_",
			  paramDomain,
			  likelihoodRoutine,
			  static_cast<void *> (&likelihoodRoutine_Data1),
			  true); // the routine computes [ln(function)]

  //------------------------------------------------------
  // SIP Step 4 of 6: Define the prior RV
  //------------------------------------------------------
  //nominal values are the ones in Forman Williams paper for detailed model
  QUESO::GslVector diagVec(paramSpace.zeroVector());
  //GAUSSIAN
  double nominal[n_params] =
  {log(3.52),log(5.06),log(1.17),log(5.75),log(1.4),-0.7,2.7,1.3,-1.4,1.0e-16,71.4,26.3,15.2,1.0e-16,249.5};
  QUESO::GslVector meanVec(paramSpace.zeroVector()); //GAUSSIAN
  for (unsigned int i = 0; i < n_params; i++) meanVec[i] = nominal[i]; //GAUSSIAN
  
  //diagVec.cwSet(.0001);
  //set prior to have st dev = 10% of nominal value
  for (unsigned int i=0; i<n_params; i++) diagVec[i] = 0.01*pow(fabs(nominal[i]),2); //GAUSSIAN
  QUESO::GslMatrix covMatrix(diagVec);
  /*
  for (unsigned int i = 0; i < 5; i++){
    covMatrix(i,i+10) = -1.0;
    covMatrix(i+10,i) = -1.0;
  }*/
  // Create a gaussian prior RV
  QUESO::GaussianVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_",paramDomain,meanVec,covMatrix);//diagVec);
   //END GAUSSIAN

  //UNIFORM
  /*
  double nominal[n_params] =
  {3.52,5.06,1.17,5.75,1.4,-0.7,2.7,1.3,-1.4,1.0e-16,71.4,26.3,15.2,1.0e-16,249.5};
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_",paramDomain);
    */
  //------------------------------------------------------
  // SIP Step 5 of 6: Instantiate the inverse problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_",  // Extra prefix before the default "rv_" prefix
           paramSpace);
        
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("",          // No extra prefix before the default "ip_" prefix
       NULL, 
       priorRv, 
       likelihoodFunctionObj, 
       postRv); 

  //------------------------------------------------------
  // SIP Step 6 of 6: Solve the inverse problem, that is,
  // set the 'pdf' and the 'realizer' of the posterior RV
  //------------------------------------------------------
  std::cout << "Solving the SIP with DRAM" 
	    << std::endl << std::endl;  
 
  //The following is set if use ip.solveWithBayesMetropolisHastings
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  for (unsigned int i = 0; i < n_params; i++) paramInitials[i] = nominal[i];

  
  //priorRv.realizer().realization(paramInitials);
  //paramInitials.cwSet(.10);

  QUESO::GslMatrix proposalCovMatrix(diagVec);
  for (unsigned int i = 0; i < n_params; i++)
    //set std dev to 10% of nominal value, same as prior
    proposalCovMatrix(i,i) = .1*diagVec[i];
  
  /*
  for (unsigned int i = 0; i < 5; i++){
    proposalCovMatrix(i,i+10) = -.0001;
    proposalCovMatrix(i+10,i) = -.0001;
  }*/
  
/* 
  for (unsigned int i = 0; i < 5; i++){
    proposalCovMatrix(i, i + 10) = 10.; 
    proposalCovMatrix(i + 10, i) = 10.; 
  }*/
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  //ip.solveWithBayesMLSampling();
  
  //================================================================
  // Statistical forward problem (SFP): find the max distance 
  // traveled by an object in projectile motion; input pdf for 'g' 
  // is the solution of the SIP above.
  //================================================================
  gettimeofday(&timevalNow, NULL);
  std::cout << "Beginning 'SFP -> Undecided QoI' at " 
            << ctime(&timevalNow.tv_sec)
            << std::endl;
	    
  //------------------------------------------------------
  // SFP Step 1 of 6: Instantiate the parameter *and* qoi spaces. 
  // SFP input RV = FIP posterior RV, so SFP parameter space
  // has been already defined.
  //------------------------------------------------------
  QUESO::VectorSpace<> qoiSpace(env, "qoi_", n_params * n_times * dim, NULL);
  //------------------------------------------------------
  // SFP Step 2 of 6: Instantiate the parameter domain 
  //------------------------------------------------------
  
  // Not necessary because input RV of the SFP = output RV of SIP. 
  // Thus, the parameter domain has been already defined.
  
  //------------------------------------------------------ 
  // SFP Step 3 of 6: Instantiate the qoi function object 
  // to be used by QUESO.
  //------------------------------------------------------
  qoiRoutine_Data qoiRoutine_Data(env, phis, tICs, times, concs, &rxnMain);
  
  QUESO::GenericVectorFunction<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    qoiFunctionObj("qoi_",
                   paramDomain,
                   qoiSpace,
                   qoiRoutine,
                   (void *) &qoiRoutine_Data);
      
  //------------------------------------------------------
  // SFP Step 4 of 6: Define the input RV
  //------------------------------------------------------
  
  // Not necessary because input RV of SFP = output RV of SIP 
  // (postRv).
      
  //------------------------------------------------------
  // SFP Step 5 of 6: Instantiate the forward problem
  //------------------------------------------------------
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix> qoiRv("qoi_", qoiSpace);
  
  QUESO::StatisticalForwardProblem<QUESO::GslVector,QUESO::GslMatrix,QUESO::GslVector,QUESO::GslMatrix>
    fp("",
       NULL,
       postRv,
       qoiFunctionObj,
       qoiRv);

  //------------------------------------------------------
  // SFP Step 6 of 6: Solve the forward problem
  //------------------------------------------------------
  std::cout << "Solving the SFP with Monte Carlo" 
            << std::endl << std::endl;  
  fp.solveWithMonteCarlo(NULL);

  //------------------------------------------------------
  gettimeofday(&timevalNow, NULL);
  if ((env.subDisplayFile()       ) && 
      (env.displayVerbosity() >= 2)) {
    *env.subDisplayFile() << "Ending run of '5 Reaction Rates with Antioch' example at "
                          << ctime(&timevalNow.tv_sec)
                          << std::endl;
  }
  if (env.fullRank() == 0) {
    std::cout << "Ending run of '5 Reaction Rates with Antioch' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  return;
}
