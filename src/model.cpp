 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This file contains the code for the forward model.
 *-----------------------------------------------------------------*/

#include "model.h"
#include "reaction_info.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <assert.h>
#include "eigen3/Eigen/Dense"
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/cea_evaluator.h>

#ifndef __EPS_ABS
#define __EPS_ABS 1e-8
#endif
#ifndef __EPS_REL
#define __EPS_REL 1e-8
#endif

// Define the function for the ODE solve
int hydrogenFunction( double t,
                      const double Y[],
                      double dYdt[],
                      void* params)
{
  // params sends the function all the reaction information
  reaction_info rxn = *(reaction_info *) params;
  //std::cout << "inside model: rxn.Reaction_set->reaction(1).forward_rate() is " << rxn.Reaction_set->reaction(0).forward_rate() << std::endl;
  unsigned int n_species = *(rxn.N_species);
  std::vector<double> molar_densities(n_species,0);
  for (unsigned int i = 0; i < n_species; i++){
    molar_densities[i] = Y[i];
    if(molar_densities[i] <= 0){
      molar_densities[i] = 0;
    }
  }

  typedef Antioch::TempCache<double> Cache;

  std::vector<double> h_RT_minus_s_R(n_species);
  rxn.Thermo->h_RT_minus_s_R(Cache(Y[n_species]),h_RT_minus_s_R);

  std::vector<double> mole_sources(n_species,0.);
  rxn.Kinetics->compute_mole_sources(
      Y[n_species], // Temperature
      molar_densities,
      h_RT_minus_s_R,
      mole_sources);

  // print chemical scheme
  /*
  const Antioch::KineticsConditions<double> cond(rxn.Temp);
  std::vector<std::vector<double>> loss_matrix;
  std::vector<std::vector<double>> prod_matrix;
  std::vector<std::vector<double>> net_matrix;

  rxn.Reaction_set->print_chemical_scheme(
      std::cout,
      cond,
      molar_densities,
      h_RT_minus_s_R,
      loss_matrix,
      prod_matrix,
      net_matrix);
      */
  
// dYdt and mole_sources are the same
  for (unsigned int i = 0; i < n_species; i++){
    if( Y[n_species] <= 1 ){
      mole_sources[i] = 0;
      std::cout << "NOTE, T dropped below 1"<<"\n\n";
    }
    dYdt[i] = mole_sources[i];
    //std::cout << "dYdt["<<i<<"] = " << dYdt[i] <<"\n";
  }

  std::vector<double> h(n_species, 0.);
  std::vector<double> cv(n_species, 0.);
  std::vector<double> cp(n_species, 0.);
  double qs = 0.;
  double Q = 0.;
  double dTdt = 0.;
  double cvtot = 0.;

  for (unsigned int s = 0; s < n_species; s++){
    h[s] = -rxn.Thermo->h(Cache(Y[n_species]),s)/rxn.Chem_mixture->R(s)* Antioch::Constants::R_universal<double>(); //negative because we want heat released
  //  std::cout<<"h = "<< h[s]<<"\n";
    //get cv and convert from mass to molar
    cv[s] = rxn.Thermo->cv(Cache(Y[n_species]),s)/rxn.Chem_mixture->R(s)* Antioch::Constants::R_universal<double>();
    cp[s] = rxn.Thermo->cp(Cache(Y[n_species]),s)/rxn.Chem_mixture->R(s)* Antioch::Constants::R_universal<double>();
    //std::cout<<"cv = "<< cv[s]<<"\n";
    qs = h[s] * mole_sources[s] * (cv[s] / cp[s]);
    Q += qs;
    cvtot += cv[s] * molar_densities[s];
  }
  dTdt = Q / cvtot;

  if( 1 && std::abs(dTdt) > 1e10 ) 
  {
    std::cout<<"Time (t) = "<<t<<"\n";
    std::cout<<"T = "<<Y[n_species]<<"\n";
    std::cout<<"She's gonna blow!!! "<<"dTdt = "<<dTdt<<"\n\n";
    /*
    for (unsigned int s = 0; s < n_species; s++){
      //std::cout << "n = " << Y[s] << std::endl;
      //std::cout << "dn_i/dt = " << mole_sources[s] << std::endl;
      //std::cout << "dq_i/dt = " << qs << std::endl;
    }*/
    //std::cout << "ntot = " << ntot << std::endl<< std::endl;
    //std::cout << "Q = " << Q << std::endl<< std::endl;
    dYdt[n_species] = dTdt;
  }
  else
    dYdt[n_species] = dTdt;
  
  //dYdt[n_species] = 0;
  //std::cout << "Y[0] = " << Y[0] << std::endl;
//  std::cout << "ntot = " << ntot << std::endl<< std::endl;
//  std::cout << "dUdt = " << Q << std::endl;
//  std::cout<<Antioch::Constants::R_universal<double>()<<std::endl;
 // std::cout<< cv*ntot<<std::endl;
 // std::cout << "dTdt = " << dTdt << std::endl<< std::endl;
  return GSL_SUCCESS;
}

//jacobian for ode solve---------------------------------------------
int hydrogenJacobian( double t, 
				const double Y[],
				double *dfdY,
				double dfdt[],
				void* params )
{
  return GSL_SUCCESS;
}

void hydrogenComputeModel(
  std::vector<double>&  initialValues,
  std::vector<double>&  timePoints,
  reaction_info*        rxn,  
  std::vector<double>&  returnValues)
{  
  // Compute model
  // GSL prep
  unsigned int dim = initialValues.size();
  unsigned int n_species = dim - 1;
  gsl_odeiv2_system sys = { hydrogenFunction, 
			   hydrogenJacobian, 
			   dim, rxn };
  
  double h = 1e-10;    // initial step-size
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45,h,1e-8,1e-4);   
  // initialize values
  double Y[dim];
  for (unsigned int i = 0; i < dim; ++i){
    Y[i] = initialValues[i];
  };
  double t = 0.0;
  double prevt;
  double sumY;
//  std::cout << "Starting Integration..." << std::endl;
  // GSL integration
  double finalTime;
  for (unsigned int i = 0; i < timePoints.size(); i++){
    finalTime = timePoints[i];
    while (t < finalTime)
      {
     //   prevt = t;
     //   sumY = 0;
        // t and y are updated and placed back into those variables
        int status = gsl_odeiv2_driver_apply( d, // necessary gsl vars
               &t,    // current time
               finalTime,    // maY time
               Y );   // current solution values (at time t)
//        std::cout<<"Time = "<<t<<"\nAfter integration Y values : \n";
//        for (unsigned int i = 0; i <= n_species; i++) std::cout<<Y[i]<<"\n";
    //    std::cout<<"t = "<<t<<"\n";
      //  std::cout<<"Y[0] = "<<Y[0]<<"\n";
        //std::cout<<"T = "<<Y[n_species]<<"\n";
//        for( i = 0; i<7;i++) sumY+=Y[i];
//        std::cout<<"N = "<<sumY<<"\n";
        // check that the evolution was successful
        #ifdef UQ_FATAL_TEST_MACRO
          UQ_FATAL_TEST_MACRO( status != GSL_SUCCESS,
             0,
             "hydrogenFiveReaction_int",
             "The status of GSL integration != GSL_SUCCESS" );
        #else 
          if ( status != GSL_SUCCESS )
          {
     //       std::cout<< "h approx = " << t-prevt<<"\n\n";
            std::cout << "ERROR: status of GSL integration != GSL_SUCCESS" <<
              std::endl;
            assert( status == GSL_SUCCESS );
          }
        #endif
      }
    //  std::cout << " h is " << h << std::endl;
    // save results, right now return values are all the species of the reduced
    // model
    for (unsigned int j = 0; j < dim; j++){
      returnValues[dim * i + j] = Y[j];}
  }
  //std::cout << "H2 = " << Y[0] << std::endl;
  //std::cout << "O2 = " << Y[1] << std::endl;
  //std::cout << "OH = " << Y[4] << std::endl;
  //std::cout << "HO2 = " << Y[5] << std::endl;
  //std::cout << "H2O = " << Y[6] << std::endl;
  //std::cout << "Temperature = " << Y[n_species] << std::endl;
  //std::cout << "conc of water is " << returnValues[6] << std::endl;
  // deallocate memory
  gsl_odeiv2_driver_free( d );
}
