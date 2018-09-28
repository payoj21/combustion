/*-------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *-------------------------------------------------------------------
 *
 */
 /*------------------------------------------------------------------
 * This file contains the code for the user defined likelihood data 
 * class and the user defined likelihood routine.
 *-----------------------------------------------------------------*/

#include "likelihood.h"
#include "reaction_info.h"
#include "model.h"
#include <cmath>
#include <stdio.h>
#include <fstream>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/kinetics_parsing.h>

// Constructor
likelihoodRoutine_Data::likelihoodRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & phis, 
    const std::vector<double> & tICs, 
    const std::vector<double> & times, 
    const std::vector<double> & concs,
    const std::vector<double> & temps,
    reaction_info * rxnInfo)
: m_env(&env),
  m_phis(phis),
  m_tICs(tICs),
  m_times(times),
  m_concs(concs),
  m_temps(temps),
  m_rxnMain(rxnInfo)
{
}

// Destructor
likelihoodRoutine_Data::~likelihoodRoutine_Data()
{
}

//------------------------------------------------------
// The user defined likelihood routine
//------------------------------------------------------

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  const QUESO::BaseEnvironment& env = *(((likelihoodRoutine_Data*) functionDataPtr)->m_env);
    
  if (paramDirection && functionDataPtr && gradVector && hessianMatrix && hessianEffect) 
  {
    // Just to eliminate INTEL compiler warnings
  }
  
  env.subComm().Barrier(); 
  
  const unsigned int n_species = 9;                 //the number of species included in the model
  const unsigned int dim = n_species + 1;           //the total dimension for the ODE solve, species + temperature
  const unsigned int n_data = 7;                    //the number of species included in the data
  const unsigned int n_times = 10;                  //the number of time points in every time series of data
  const unsigned int n_phis = 1;                  //the number of time points in every time series of data
  const unsigned int n_tICs = 1;                  //the number of time points in every time series of data
  const unsigned int n_reactions = 5;               //the number of reactions in the reduced model
  const unsigned int n_params = 3 * n_reactions;    //the number of parameters in the reduced model (A, n, E), k = A T^n exp(-E/RT)

  // Compute likelihood 
  // get data from likelihood data structure
  const std::vector<double>&  phis	
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_phis;
  const std::vector<double>& 	times
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_times;
  const std::vector<double>& 	tICs
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_tICs;
  const std::vector<double>& 	concs
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_concs;
  const std::vector<double>&  temps
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_temps;
  reaction_info *        rxn5
    = ((likelihoodRoutine_Data *) functionDataPtr)->m_rxnMain;

  double variance = .01;
  //double variance = .001; // what the variance was usually set to

  //set up lambda vector for loop, right now just one
  std::vector<double> phiPoints(n_phis,0.);
  for (unsigned int j = 0; j < n_phis; j++){
    phiPoints[j] = phis[n_tICs * n_times * j];
  }
  std::vector<double> tempPoints(n_tICs,0.);
  for (unsigned int j = 0; j < n_tICs; j++){
    tempPoints[j] = tICs[n_times * j];
  }

  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i];
  }

  //*****define initial conditions = state + temp
  std::vector<double> initial_conditions(dim, 0.0);
  //get initial densities
  //H_2 will be set within upcoming loop over lambda
  initial_conditions[1] = 1.0; //O2
  initial_conditions[8] = 3.78; //N2

  //return 10 time points of H2O
  std::vector<double> returnValues(n_times * dim, 0.);

  //temporary values to reset rates with antioch
  std::vector<double> tempValues(n_params);
  //GAUSSIAN
   //transform gaussian to lognormal and scale
  for (unsigned int i = 0; i < n_reactions; i++) tempValues[i] = exp(paramValues[i]) * rxn5->Scales[i];
  //multiply newValues for A by constants to keep params close to 1
  for (unsigned int i = n_reactions; i < n_params; i++) tempValues[i] = paramValues[i] * rxn5->Scales[i];
  //END GAUSSIAN

  //std::cout << "rxn5.Temp is " << rxn5.Temp << std::endl;
  //std::cout << "rxn5.Reaction_set->reaction(1) is " << rxn5.Reaction_set->reaction(1) << std::endl;
  //int n = 0;
  //std::cout << "before: rxn5->Reaction_set->reaction(n).forward_rate() is " << rxn5->Reaction_set->reaction(n).forward_rate() << std::endl;

  std::vector<double> newValues(3);
  for (unsigned int i = 0; i < n_reactions; i++) {
    newValues[0] = tempValues[i];
    newValues[1] = tempValues[i + n_reactions];
    newValues[2] = tempValues[i + 2*n_reactions];
    reset_rate(rxn5->Reaction_set->reaction(i).forward_rate(),newValues);
  }

  //std::cout << " after: rxn5->Reaction_set->reaction(n).forward_rate() is " << rxn5->Reaction_set->reaction(n).forward_rate() << std::endl;
  /*
  std::cout<< "after\n";
  rxn5->Reaction_set->print_chemical_scheme(
      std::cout,
      cond,
      molar_densities,
      rxn5->H_RT_minus_s_R,
      loss_matrix,
      prod_matrix,
      net_matrix);
*/

  double misfitValue = 0.;
  double diff = 0.;

  for (unsigned int i = 0; i < n_phis; i++) {
   for (unsigned int ii = 0; ii < n_tICs; ii++) {
    initial_conditions[0] = 2. * phiPoints[i]; //amount of H_2
    initial_conditions[n_species] = tempPoints[ii];
    try
     {
      hydrogenComputeModel(initial_conditions,timePoints,rxn5,returnValues);
//      std::cout << "Finished compute model" << std::endl;
      for (unsigned int j = 0; j < n_times; j++){
        for (unsigned int k = 0; k < n_data; k++){ //H2O2 & N2 not included in data
          diff = (returnValues[dim * j + k] - concs[(n_tICs*n_times*n_species) * i + (n_times*n_species) *ii + n_species*j + k]);
          //std::cout<<"likelihood: " << (n_tICs*n_times*n_species) * i + (n_times*n_species) *ii + n_species*j + k << "\n";
          misfitValue += diff * diff / variance;
        }
        diff = returnValues[dim * (j+1) - 1] - temps[(n_times*n_tICs)*i + n_times*ii + j];
        //std::cout << returnValues[dim * (j+1) - 1] << " and " << temps[j] << "\n";
        misfitValue += diff * diff / 1000; //assume variance on temp data is 1000
      }
     } catch( int exception )
     {
      misfitValue = 1000000;
   }

    }
   } 

//  std::cout << " the misfit is " << misfitValue << std::endl;
  return (-0.5 * misfitValue);  
}
