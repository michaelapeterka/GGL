#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include "probe_energiecalculation.h"

using namespace std;

int main()
{
vector<string> probe_smiles;

probe_smiles.push_back("[H:1][C:1]#[N:1]");

double energie_producedSmiles = energie_calculation(probe_smiles);

cout << "Energy_producedSmiles" << energie_producedSmiles << endl; 
return 0;
}

