#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include "probe_energiecalculation.h"
#include <map>

using namespace std;

int main()
{
  //Probeversion/Probewerte etc.
  std::map<std::string,double> probe_smiles;
  std::map<std::string,double> metabolites;
  std::string smiles = "[H:1][C:1]#[N:1]";
  std::string m_smiles = "[H:1]";
  std::string m_smiles1 = "[C:1][N:1]";

  metabolites[m_smiles]= 54.0;
  probe_smiles[smiles] = 5.6978;
  metabolites[m_smiles1] = 4565.09; 
  
  
double energie_producedSmiles = energie_calculation(probe_smiles);
double energie_metabolites = energie_calculation(metabolites);
 double a0 = energie_producedSmiles - energie_metabolites;
cout << "Energy_producedSmiles:" << energie_producedSmiles << endl;
 cout << "Energy_metabolites: " << energie_metabolites << endl;
 cout <<"a0 = " << a0 << endl;

 
return 0;
}

