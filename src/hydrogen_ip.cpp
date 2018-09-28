 /*------------------------------------------------------------------
 * Brief description of this file: 
 * 
 * This is the inverse problem to calibrate the 5 rate constants for the 5
 * reaction reduced model of hydrogen combustion.
 *
 * Right now we don't have a set QoI, but there will probably be one in the near
 * future.
 * 
 * The code consists of 7 files:
 * - 'hydrogen_5.C' (this file)
 * - 'compute.C' (the driving application code)
 * - 'compute.h'
 * - 'likelihood.C' (necessary for the SIP)
 * - 'likelihood.h'
 * - 'qoi.C' (necessary for the SFP)
 * - 'qoi.h'
 *-----------------------------------------------------------------*/

//#include <queso/Environment.h>
#include <compute.h>
int main(int argc, char* argv[])
{
  // Initialize QUESO environment

#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(argv[1],"",NULL);
#endif
  // Call application
  computeReactionRates(*env);

  // Finalize QUESO environment
  delete env;

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif

  return 0;
}
