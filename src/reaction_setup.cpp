#include "reaction_setup.h"

class hydrogen_setup
{
    public:
        hydrogen_setup(std::string data_input_name);
        ~hydrogen_setup();
    
        
        std::vector<std::string> species{
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
            return species_str_list;
        }
        std::vector<double> data_phis{

};

hydrogen_setup::hydrogen_setup(std::string data_input_name)
{

  FILE *dataFile;
  dataFile = fopen(data_input_name,"r");

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

  std::vector<double> data_phis;
  std::vector<double> data_initial_temps;
  std::vector<double> data_times;
  std::vector<double> data_concs;
  std::vector<double> data_temps;
  std::cout << "HELLO? " << std::endl;
  while (fscanf(dataFile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tmpPhis, &tmpTICs,
        &tmpTimes, &tmpH2, &tmpO2, &tmpH, &tmpO, &tmpOH, &tmpHO2, &tmpH2O, &tmpH2O2, &tmpN2, &tmpTemps) != EOF) {
    std::cout << "HELLO : " << numLines << std::endl;
    data_phis.push_back(tmpPhis);
    data_initial_temps.push_back(tmpTICs);
    data_times.push_back(tmpTimes);
    data_concs.push_back(tmpH2);
    data_concs.push_back(tmpO2);
    data_concs.push_back(tmpH);
    data_concs.push_back(tmpO);
    data_concs.push_back(tmpOH);
    data_concs.push_back(tmpHO2);
    data_concs.push_back(tmpH2O);
    data_concs.push_back(tmpH2O2);
    data_concs.push_back(tmpN2);
    data_temps.push_back(tmpTemps);
    numLines++;
    numData+=n_species;
  }

}
