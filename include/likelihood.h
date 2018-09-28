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
 * Brief description of this file: 
 *
 * This is the header file for hydrogen_likelihood.C. 
 *-----------------------------------------------------------------*/

#ifndef __HYDROGEN_LIKELIHOOD_H__
#define __HYDROGEN_LIKELIHOOD_H__

#include "reaction_info.h"
#include <queso/GslMatrix.h>

struct likelihoodRoutine_Data // user defined class
{
  likelihoodRoutine_Data(
      const QUESO::BaseEnvironment& env,
      const std::vector<double> & phis,
      const std::vector<double> & tICs,
      const std::vector<double> & times, 
      const std::vector<double> & concs, 
      const std::vector<double> & temps, 
      reaction_info * rxnInfo);
 ~likelihoodRoutine_Data();

  const QUESO::BaseEnvironment* m_env;
  const std::vector<double> & m_phis;   //
  const std::vector<double> & m_tICs;   //
  const std::vector<double> & m_times;   //
  const std::vector<double> & m_concs;   // 
  const std::vector<double> & m_temps;   //
  reaction_info       * m_rxnMain;
};

double likelihoodRoutine( // user defined routine
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect);

#endif
