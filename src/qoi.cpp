/*------------------------------------------------------------------------
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
 *------------------------------------------------------------------------
 *
 */
 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the user defined qoi routine.
 *-----------------------------------------------------------------*/

#include "qoi.h"
#include "model.h"
#include "reaction_info.h"
#include <cmath>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/kinetics_parsing.h>
//------------------------------------------------------
/// The actual (user-defined) qoi routine
//------------------------------------------------------
// Constructor
qoiRoutine_Data::qoiRoutine_Data(
    const QUESO::BaseEnvironment& env,
    const std::vector<double> & phis,
    const std::vector<double> & tICs,
    const std::vector<double> & times,
    std::vector<double> & concs,
    reaction_info * rxnInfo)
: m_env(&env),
  m_phis(phis),
  m_tICs(tICs),
  m_times(times),
  m_concs(concs),
  m_rxnMain(rxnInfo)
{
}

// Destructor
qoiRoutine_Data::~qoiRoutine_Data()
{
}

void
qoiRoutine(
  const QUESO::GslVector&                    paramValues,
  const QUESO::GslVector*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslVector&                    qoiValues,
        QUESO::DistArray<QUESO::GslVector*>* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  const QUESO::BaseEnvironment& env = paramValues.env();

  if (paramDirection && 
      gradVectors    &&
      hessianEffects &&
      hessianMatrices) {
    // Logic just to avoid warnings from INTEL compiler
  }

  const unsigned int n_species = 9;                 //the number of species included in the model
  const unsigned int dim = n_species + 1;           //the total dimension for the ODE solve, species + temperature
  const unsigned int n_data = 7;                    //the number of species included in the data
  const unsigned int n_times = 10;                  //the number of time points in every time series of data
  const unsigned int n_phis = 1;                    //the number of initial conditions (phi = equivalence ratio) 
  const unsigned int n_tICs = 1;                    //the number of initial conditions (phi = equivalence ratio) 
  const unsigned int n_reactions = 5;               //the number of reactions in the reduced model
  const unsigned int n_params = 3 * n_reactions;    //the number of parameters in the reduced model (A, n, E), k = A T^n exp(-E/RT)
  
  // get data from qoi data structure
  const std::vector<double>&  phis	
    = ((qoiRoutine_Data *) functionDataPtr)->m_phis;
  const std::vector<double>&  tICs	
    = ((qoiRoutine_Data *) functionDataPtr)->m_tICs;
  const std::vector<double>& 	times
    = ((qoiRoutine_Data *) functionDataPtr)->m_times;
  const std::vector<double>& 	concs
    = ((qoiRoutine_Data *) functionDataPtr)->m_concs;
  reaction_info *       rxn5
    = ((qoiRoutine_Data *) functionDataPtr)->m_rxnMain;
  
  double variance = .001;

  //define initial conditions = state + temp
  std::vector<double> initial_conditions(dim, 0.0);
  //get initial densities

  //H_2 will be set within upcoming loop over lambda
  initial_conditions[1] = 1.0; //O2
  initial_conditions[8] = 3.78; //N2

  //set up lambda vector for loop, right now just one
  std::vector<double> phiPoints(n_phis,0.);
  for (unsigned int j = 0; j < n_phis; j++){
    phiPoints[j] = phis[n_tICs * n_times * j];
  }
  std::vector<double> tempPoints(n_tICs,0.);
  for (unsigned int i = 0; i < n_tICs; i++){
    tempPoints[i] = tICs[n_times * i];
  }
  std::vector<double> timePoints(n_times,0.);
  for (unsigned int i = 0; i < n_times; i++){
    timePoints[i] = times[i];
  }

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

  //UNIFORM
  /*
  //multiply newValues for A by constants to keep params close to 1
  for (unsigned int i = 0; i < n_params; i++) tempValues[i] = paramValues[i] * scaling[i];
  */ //END UNIFORM

  //reset the rates in kinetics model
  std::vector<double> newValues(3);
  for (unsigned int i = 0; i < n_reactions; i++) {
    newValues[0] = tempValues[i];
    newValues[1] = tempValues[i + n_reactions];
    newValues[2] = tempValues[i + 2*n_reactions];
    reset_rate(rxn5->Reaction_set->reaction(i).forward_rate(),newValues);
  }

  for (unsigned int i = 0; i < n_phis; ++i) {
   for (unsigned int ii = 0; ii < n_tICs; ii++){
    initial_conditions[0] = 2. * phiPoints[i];  //amount of H_2
    initial_conditions[n_species] = tempPoints[ii];
    hydrogenComputeModel(initial_conditions,timePoints,rxn5,returnValues);
    //std::cout<< "qoi: ret val = " <<  returnValues[7 * 9 + 6] << std::endl; 
    for (unsigned int j = 0; j < returnValues.size(); j++){
        //std::cout << "j = " << j << " and k = " << k << "\n\n";
        qoiValues[(n_tICs*n_times*dim) * i + (n_times*dim) *ii + j] = returnValues[j];
      }
    }
  }

// for (unsigned int i = 0; i < 90; i++) qoiValues[i] = i;
  return;
}
